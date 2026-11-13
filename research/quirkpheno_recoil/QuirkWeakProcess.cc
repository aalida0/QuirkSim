#include "quirkcore/QuirkWeakProcess.h"
#include "quirkcore/QuirkStringStore.h"
#include "quirkcore/QuirkTrackStatus.h"
#include "quirkcore/TwoBodyKinematics.h"
#include "quirkcore/DirectionControl.h"
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
  if (!(QPWeakConfig::tauPs(m_neutral->GetPDGEncoding())>0) || (QPWeakConfig::directionSign()!=1 && QPWeakConfig::directionSign()!=-1))
    G4Exception("QuirkWeakProcess", "BadWeakPilot", FatalException, "positive threshold and +/-1 direction required");
}
G4double QuirkWeakProcess::PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition* condition) {
  *condition=StronglyForced;
  return DBL_MAX; // transportation caps the localized clock; legacy mode remains post-step
}
G4VParticleChange* QuirkWeakProcess::PostStepDoIt(const G4Track& track, const G4Step&) {
  m_change.Initialize(track);
  if (track.GetTrackStatus()==fStopAndKill || track.GetTrackStatus()==fKillTrackAndSecondaries ||
      track.GetProperTime() + (QPWeakConfig::localizedClock() ? 1e-9*CLHEP::picosecond : 0.) < QPWeakConfig::tauPs(track.GetDefinition()->GetPDGEncoding())*CLHEP::picosecond) return &m_change;
  auto& store=QuirkStringStore::instance();
  auto& string=store.stringFor(track.GetDefinition()->GetPDGEncoding());
  if (!string.IsSourceAlive() || !string.IsSourceInitialized())
    G4Exception("QuirkWeakProcess", "WeakOrphan", FatalException, "decay without a live initialized partner");
  if (QPWeakConfig::generalPlan()) return generalDecay(track);
  const int pionPDG=track.GetDefinition()->GetPDGCharge()>0 ? 211 : -211;
  auto* pion=G4ParticleTable::GetParticleTable()->FindParticle(pionPDG);
  if (!pion) G4Exception("QuirkWeakProcess", "MissingPion", FatalException, "pion not registered");
  const auto p=track.GetMomentum();
  const auto direction=QPWeak::controlledDirection(QPWeakConfig::directionMode(),QPWeakConfig::collinearSign(),
      {p.x(),p.y(),p.z()},QPWeakConfig::pionDirection(track.GetDefinition()->GetPDGEncoding()),
      track.GetDynamicParticle()->GetMass(),QPWeakConfig::controlBoostZ());
  auto k=QPWeak::pionDecay(track.GetDynamicParticle()->GetMass(), m_neutral->GetPDGMass(),
                          pion->GetPDGMass(), {p.x(),p.y(),p.z()}, direction);
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
    << ",\"direction_mode\":\"" << QPWeakConfig::directionMode() << "\""
    << ",\"control_boost_z\":" << QPWeakConfig::controlBoostZ()
    << ",\"rest_frame_direction\":[" << direction[0] << ',' << direction[1] << ',' << direction[2] << ']'
    << ",\"zero_momentum_axis_fallback\":" << (QPWeakConfig::directionMode()=="collinear" && p.mag2()==0 ? "true" : "false")
    << ",\"time_ps\":" << track.GetGlobalTime()/CLHEP::picosecond
    << ",\"position_mm\":[" << track.GetPosition().x()/CLHEP::mm << "," << track.GetPosition().y()/CLHEP::mm << "," << track.GetPosition().z()/CLHEP::mm << "]"
    << ",\"proper_ps\":" << track.GetProperTime()/CLHEP::picosecond
    << ",\"clock_mode\":\"" << (QPWeakConfig::localizedClock() ? "localized" : "legacy") << "\""
    << ",\"threshold_ps\":" << QPWeakConfig::tauPs(track.GetDefinition()->GetPDGEncoding()) << ",\"residual_mev\":" << residual;
  auto vector=[&](const char* name,const QPWeak::Four& v) {
    s << ",\"" << name << "\":[" << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << ']';
  };
  vector("parent_mev",k.parent); vector("daughter_mev",k.daughter); vector("pion_mev",k.pion);
  s << '}'; G4cout << s.str() << G4endl;
  return &m_change;
}

G4VParticleChange* QuirkWeakProcess::generalDecay(const G4Track& track) {
  const int pdg=track.GetDefinition()->GetPDGEncoding();
  const auto& plan=QPWeakConfig::generalPlan()->endpoint[pdg>0?0:1];
  const auto p=track.GetMomentum();const double M=track.GetDynamicParticle()->GetMass();
  if(std::abs(M-plan.M)>1e-6*CLHEP::MeV || std::abs(m_neutral->GetPDGMass()-plan.m)>1e-6*CLHEP::MeV)
    G4Exception("QuirkWeakProcess","WeakPlanMass",FatalException,"inclusive plan and endpoint masses disagree");
  const QPWeak::Four parent{std::sqrt(M*M+p.mag2()),p.x(),p.y(),p.z()};
  std::array<double,3> axis{};
  for(size_t j=1;j<plan.products.size();++j)for(int k=0;k<3;++k)axis[k]+=plan.products[j].p4[k+1];
  const auto target=QPWeak::controlledDirection(QPWeakConfig::directionMode(),QPWeakConfig::collinearSign(),
      {p.x(),p.y(),p.z()},axis,M,QPWeakConfig::controlBoostZ());
  auto products=QPWeak::labProducts(plan,parent,QPWeakConfig::directionMode()=="planned"?nullptr:&target);
  QPWeak::Four total{};double residual=0.;
  for(const auto& v:products)for(int k=0;k<4;++k)total[k]+=v.p4[k];
  for(int k=0;k<4;++k)residual=std::max(residual,std::abs(total[k]-parent[k]));
  if(residual>1e-6*CLHEP::MeV)
    G4Exception("QuirkWeakProcess","WeakGeneralMomentum",FatalException,"inclusive weak momentum residual exceeds tolerance");
  std::vector<G4Track*> children;
  for(size_t j=0;j<products.size();++j) {
    const auto& v=products[j];const auto* def=j==0?m_neutral:G4ParticleTable::GetParticleTable()->FindParticle(v.pdg);
    if(!def)G4Exception("QuirkWeakProcess","WeakProductMissing",FatalException,"inclusive weak product is not registered");
    // Pin the explicitly recorded mass convention. Geant4 particle tables can
    // differ from the plan's PDG constants by fractions of a keV.
    if(std::abs(v.mass-def->GetPDGMass())>1e-3*CLHEP::MeV)
      G4Exception("QuirkWeakProcess","WeakProductMass",FatalException,"plan/table mass difference exceeds 1 keV");
    auto* dynamic=new G4DynamicParticle(def,G4ThreeVector(v.p4[1],v.p4[2],v.p4[3]));
    dynamic->SetMass(v.mass);dynamic->SetKineticEnergy(v.p4[0]-v.mass);
    const auto actual=dynamic->GetMomentum();
    if(std::max({std::abs(actual.x()-v.p4[1]),std::abs(actual.y()-v.p4[2]),std::abs(actual.z()-v.p4[3])})>1e-6*CLHEP::MeV)
      G4Exception("QuirkWeakProcess","WeakDynamicMomentum",FatalException,"Geant4 dynamic particle differs from planned four-vector");
    auto* child=new G4Track(dynamic,track.GetGlobalTime(),track.GetPosition());
    child->SetTouchableHandle(track.GetTouchableHandle());
    if(j==0)child->SetProperTime(track.GetProperTime());
    children.push_back(child);
  }
  auto& store=QuirkStringStore::instance();
  store.beginHandoff(pdg,track.GetTrackID(),m_neutral->GetPDGEncoding());
  m_change.SetNumberOfSecondaries(static_cast<G4int>(children.size()));
  for(auto* child:children)m_change.AddSecondary(child);
  m_change.ProposeEnergy(0.);m_change.ProposeTrackStatus(fStopAndKill);
  store.stringFor(pdg).SetKineticEnergy(products[0].p4[0]-m_neutral->GetPDGMass());
  std::ostringstream s;
  s<<std::setprecision(17)<<"QP_WEAK_JSON {\"kind\":\"decay\",\"schema\":2,\"channel\":\""<<plan.channel<<"\",\"parent_pdg\":"<<pdg
    <<",\"parent_track\":"<<track.GetTrackID()<<",\"daughter_pdg\":"<<m_neutral->GetPDGEncoding()
    <<",\"parent_mass_mev\":"<<M<<",\"daughter_mass_mev\":"<<m_neutral->GetPDGMass()
    <<",\"direction_mode\":\""<<QPWeakConfig::directionMode()<<"\",\"control_boost_z\":"<<QPWeakConfig::controlBoostZ()
    <<",\"time_ps\":"<<track.GetGlobalTime()/CLHEP::picosecond<<",\"proper_ps\":"<<track.GetProperTime()/CLHEP::picosecond
    <<",\"threshold_ps\":"<<plan.properPs<<",\"residual_mev\":"<<residual
    <<",\"position_mm\":["<<track.GetPosition().x()/CLHEP::mm<<','<<track.GetPosition().y()/CLHEP::mm<<','<<track.GetPosition().z()/CLHEP::mm<<']';
  auto vec=[&](const char* name,const QPWeak::Four& v){s<<",\""<<name<<"\":["<<v[0]<<','<<v[1]<<','<<v[2]<<','<<v[3]<<']';};
  vec("parent_mev",parent);vec("daughter_mev",products[0].p4);s<<",\"products\":[";
  for(size_t j=1;j<products.size();++j) {
    if(j>1)s<<',';
    const auto& v=products[j];s<<"{\"pdg\":"<<v.pdg<<",\"mass_mev\":"<<v.mass<<",\"p4_mev\":["<<v.p4[0]<<','<<v.p4[1]<<','<<v.p4[2]<<','<<v.p4[3]<<"]}";
  }
  s<<"]}";G4cout<<s.str()<<G4endl;
  return &m_change;
}
