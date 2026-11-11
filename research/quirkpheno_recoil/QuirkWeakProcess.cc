#include "quirkcore/QuirkWeakProcess.h"
#include "quirkcore/QuirkStringStore.h"
#include "quirkcore/QuirkTrackStatus.h"
#include "quirkcore/TwoBodyKinematics.h"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4Exception.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include <algorithm>
#include <iomanip>
#include <sstream>

QuirkWeakProcess::QuirkWeakProcess(const G4ParticleDefinition* neutral)
 : G4VDiscreteProcess("QuirkWeakPilot", fDecay), m_neutral(neutral) {
  if (!(QPWeakConfig::tauPs()>0) || (QPWeakConfig::directionSign()!=1 && QPWeakConfig::directionSign()!=-1))
    G4Exception("QuirkWeakProcess", "BadWeakPilot", FatalException, "positive threshold and +/-1 direction required");
}
G4double QuirkWeakProcess::PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition* condition) {
  *condition=StronglyForced;
  return DBL_MAX; // transportation caps the localized clock; legacy mode remains post-step
}
G4VParticleChange* QuirkWeakProcess::PostStepDoIt(const G4Track& track, const G4Step&) {
  m_change.Initialize(track);
  if (track.GetTrackStatus()==fStopAndKill || track.GetTrackStatus()==fKillTrackAndSecondaries ||
      track.GetProperTime() + (QPWeakConfig::localizedClock() ? 1e-9*CLHEP::picosecond : 0.) < QPWeakConfig::tauPs()*CLHEP::picosecond) return &m_change;
  auto& store=QuirkStringStore::instance();
  auto& string=store.stringFor(track.GetDefinition()->GetPDGEncoding());
  if (!string.IsSourceAlive() || !string.IsSourceInitialized())
    G4Exception("QuirkWeakProcess", "WeakOrphan", FatalException, "decay without a live initialized partner");
  const int pionPDG=track.GetDefinition()->GetPDGCharge()>0 ? 211 : -211;
  auto* pion=G4ParticleTable::GetParticleTable()->FindParticle(pionPDG);
  if (!pion) G4Exception("QuirkWeakProcess", "MissingPion", FatalException, "pion not registered");
  const auto p=track.GetMomentum();
  auto k=QPWeak::pionDecay(track.GetDynamicParticle()->GetMass(), m_neutral->GetPDGMass(),
                          pion->GetPDGMass(), {p.x(),p.y(),p.z()}, {0,QPWeakConfig::directionSign(),0});
  double residual=0;
  for (int i=0;i<4;++i) residual=std::max(residual,std::abs(k.parent[i]-k.daughter[i]-k.pion[i]));
  if (residual>1e-6*CLHEP::MeV)
    G4Exception("QuirkWeakProcess", "WeakMomentum", FatalException, "local four-momentum residual exceeds 1e-6 MeV");
  auto make=[&](const G4ParticleDefinition* def, const QPWeak::Four& v) {
    auto* dp=new G4DynamicParticle(def,G4ThreeVector(v[1],v[2],v[3]));
    auto* child=new G4Track(dp,track.GetGlobalTime(),track.GetPosition());
    child->SetTouchableHandle(track.GetTouchableHandle());
    return child;
  };
  auto* daughter=make(m_neutral,k.daughter);
  // This is a continuing endpoint worldline, so retain its accumulated proper clock.
  daughter->SetProperTime(track.GetProperTime());
  auto* emitted=make(pion,k.pion);
  store.beginHandoff(track.GetDefinition()->GetPDGEncoding(),track.GetTrackID(),m_neutral->GetPDGEncoding());
  m_change.SetNumberOfSecondaries(2);
  m_change.AddSecondary(daughter);
  m_change.AddSecondary(emitted);
  m_change.ProposeEnergy(0.);
  m_change.ProposeTrackStatus(fStopAndKill);
  string.SetKineticEnergy(k.daughter[0]-m_neutral->GetPDGMass());
  std::ostringstream s;
  s << std::setprecision(17) << "QP_WEAK_JSON {\"kind\":\"decay\",\"parent_track\":" << track.GetTrackID()
    << ",\"parent_pdg\":" << track.GetDefinition()->GetPDGEncoding() << ",\"daughter_pdg\":" << m_neutral->GetPDGEncoding()
    << ",\"pion_pdg\":" << pionPDG
    << ",\"parent_mass_mev\":" << track.GetDynamicParticle()->GetMass()/CLHEP::MeV
    << ",\"daughter_mass_mev\":" << m_neutral->GetPDGMass()/CLHEP::MeV
    << ",\"pion_mass_mev\":" << pion->GetPDGMass()/CLHEP::MeV
    << ",\"time_ps\":" << track.GetGlobalTime()/CLHEP::picosecond
    << ",\"position_mm\":[" << track.GetPosition().x()/CLHEP::mm << "," << track.GetPosition().y()/CLHEP::mm << "," << track.GetPosition().z()/CLHEP::mm << "]"
    << ",\"proper_ps\":" << track.GetProperTime()/CLHEP::picosecond
    << ",\"clock_mode\":\"" << (QPWeakConfig::localizedClock() ? "localized" : "legacy") << "\""
    << ",\"threshold_ps\":" << QPWeakConfig::tauPs() << ",\"residual_mev\":" << residual;
  auto vector=[&](const char* name,const QPWeak::Four& v) {
    s << ",\"" << name << "\":[" << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << ']';
  };
  vector("parent_mev",k.parent); vector("daughter_mev",k.daughter); vector("pion_mev",k.pion);
  s << '}'; G4cout << s.str() << G4endl;
  return &m_change;
}
