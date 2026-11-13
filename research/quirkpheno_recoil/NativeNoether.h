#ifndef QP_NATIVE_NOETHER_H
#define QP_NATIVE_NOETHER_H
// Opt-in, single-event vacuum/uniform-Bz diagnostic. No transport modifications.
#include "quirkcore/NoetherCharges.h"
#include "quirkcore/EncounterRecorder.h"
#include "quirkcore/QuirkStringStore.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace QPNativeNoether {
using namespace QPNoether;
inline bool enabled(){const auto* p=std::getenv("QP_NOETHER_OUTPUT");return p && *p;}
inline V position(const G4ThreeVector& x,double time){return {CLHEP::c_light*time/CLHEP::mm,x.x()/CLHEP::mm,x.y()/CLHEP::mm,x.z()/CLHEP::mm};}
inline V momentum(const G4LorentzVector& p){return {p.e()/CLHEP::MeV,p.px()/CLHEP::MeV,p.py()/CLHEP::MeV,p.pz()/CLHEP::MeV};}
struct State {
  std::array<V,2> x{},p{},initialX{},initialP{};
  std::array<bool,2> seen{};std::array<int,2> skipTrack{};
  std::array<double,2> properPs{},mass{};std::array<int,2> pdg{};
  QPEncounter::Recorder encounters;std::ofstream encounterOutput;
  Charge flux{};long double torqueZ=0;unsigned long steps=0,snapshots=0,vertices=0;
  std::ofstream output;int event=-1;
};
inline State& state(){static thread_local State s;return s;}
inline void event(){
  auto& s=state();const auto* e=G4EventManager::GetEventManager()->GetConstCurrentEvent();
  if(!e)throw std::runtime_error("Noether: no active event");
  if(s.event==e->GetEventID())return;
  if(s.event!=-1)throw std::runtime_error("Noether diagnostic requires one event per process");
  s.event=e->GetEventID();s.output.open(std::getenv("QP_NOETHER_OUTPUT"));s.output<<std::setprecision(21);
  if(!s.output)throw std::runtime_error("Noether output could not be opened");
}
inline void move(int i,const V& next,long double charge){
  auto& s=state();const char* field=std::getenv("QP_NOETHER_BZ_T");
  if(!field)throw std::runtime_error("Noether requires declared uniform Bz");
  const long double B=std::stold(field);
  // q c B in MeV/mm. Integrates Lorentz axial torque exactly in uniform Bz.
  auto r2=[](const V& x){return x[1]*x[1]+x[2]*x[2];};
  s.torqueZ-=.299792458L*B*charge*(r2(next)-r2(s.x[i]))/2;
  s.x[i]=next;
}
inline void snapshot(const char* reason){
  auto& s=state();auto& store=QuirkStringStore::instance();
  if(!s.seen[0] || !s.seen[1])return;
  std::array<std::vector<V>,2> rays;
  long double tension=0;
  for(int i=0;i<2;++i){auto& f=store.stringFor(i==0?17:-17);
    if(!f.NoetherReady())return;
    tension=f.GetStringForce()/(CLHEP::MeV/CLHEP::mm);
    for(const auto& r:f.GetStringVectors())rays[i].push_back({r.t()/CLHEP::mm,r.x()/CLHEP::mm,r.y()/CLHEP::mm,r.z()/CLHEP::mm});
  }
  auto q=total(s.x,s.p,rays,tension);auto close=closure(s.x,rays);
  V initial=add(s.initialP[0],s.initialP[1]);T initialJ=add(wedge(s.initialX[0],s.initialP[0]),wedge(s.initialX[1],s.initialP[1]));
  auto& o=s.output;auto vector=[&](const auto& a){o<<'[';for(size_t k=0;k<a.size();++k)o<<(k?",":"")<<a[k];o<<']';};
  o<<"{\"schema\":1,\"reason\":\""<<reason<<"\",\"step\":"<<s.steps<<",\"vertices\":"<<s.vertices<<",\"endpoint_ct_mm\":["<<s.x[0][0]<<','<<s.x[1][0]<<"],\"contour_closure_mm\":";vector(close);
  o<<",\"mechanical_p_mev\":";vector(q.p);o<<",\"mechanical_j_mev_mm\":";vector(q.j);
  o<<",\"emitted_p_mev\":";vector(s.flux.p);o<<",\"emitted_j_mev_mm\":";vector(s.flux.j);
  o<<",\"initial_p_mev\":";vector(initial);o<<",\"initial_j_mev_mm\":";vector(initialJ);
  o<<",\"endpoint_proper_ps\":";vector(s.properPs);o<<",\"endpoint_pdg\":";vector(s.pdg);
  o<<",\"external_torque_z_mev_mm\":"<<s.torqueZ<<",\"string_p_mev\":";vector(add(q.p,scale(add(s.p[0],s.p[1]),-1)));
  o<<"}\n";o.flush();if(!o)throw std::runtime_error("Noether output write failed");++s.snapshots;
}
inline void observe(int i,const V& x,const V& p,long double mass,int pdg){
  const char* path=std::getenv("QP_ENCOUNTER_OUTPUT");if(!path || !*path)return;
  auto& s=state();
  if(!s.encounterOutput.is_open()){s.encounterOutput.open(path);s.encounterOutput<<std::setprecision(21);if(!s.encounterOutput)throw std::runtime_error("encounter output open failed");}
  const char* btext=std::getenv("QP_WEAK_CONTROL_BOOST_Z");if(!btext)throw std::runtime_error("encounters need declared CM boost");
  const long double beta=std::stold(btext),gamma=1/std::sqrt(1-beta*beta);
  QPEncounter::Point a;a.t=gamma*(x[0]-beta*x[3])/.299792458L;
  a.x={x[1],x[2],gamma*(x[3]-beta*x[0])};a.p={p[1],p[2],gamma*(p[3]-beta*p[0])};a.e=gamma*(p[0]-beta*p[3]);a.mass=mass;a.pdg=pdg;
  s.encounters.push(i,a,[&](const QPEncounter::Encounter& e){
    auto& o=s.encounterOutput;o<<"{\"kind\":\"return\",\"time_cm_ps\":"<<e.t<<",\"separation_mm\":"<<e.miss<<",\"endpoints\":[";
    for(int j=0;j<2;++j){if(j)o<<',';const auto& q=e.endpoint[j];o<<"{\"pdg\":"<<q.pdg<<",\"mass_mev\":"<<q.mass<<",\"x_mm\":["<<q.x[0]<<','<<q.x[1]<<','<<q.x[2]<<"],\"p4_mev\":["<<q.e<<','<<q.p[0]<<','<<q.p[1]<<','<<q.p[2]<<"]}";}
    o<<"]}\n";if(!o)throw std::runtime_error("encounter output write failed");
  });
}
inline void step(const G4Step* step){
  if(!enabled())return;
  event();auto& s=state();const auto* t=step->GetTrack();
  const int i=t->GetDefinition()->GetPDGEncoding()>0?0:1;
  const auto* pre=step->GetPreStepPoint();
  if(!s.seen[i]){s.seen[i]=true;s.x[i]=s.initialX[i]=position(pre->GetPosition(),pre->GetGlobalTime());
    s.initialP[i]=momentum(G4LorentzVector(pre->GetMomentum(),pre->GetTotalEnergy()));
    observe(i,s.initialX[i],s.initialP[i],t->GetDynamicParticle()->GetMass()/CLHEP::MeV,t->GetDefinition()->GetPDGEncoding());}
  if(s.skipTrack[i]==t->GetTrackID()){s.skipTrack[i]=0;return;}
  move(i,position(t->GetPosition(),t->GetGlobalTime()),t->GetDynamicParticle()->GetCharge()/CLHEP::eplus);
  s.p[i]=momentum(t->GetDynamicParticle()->Get4Momentum());
  s.properPs[i]=t->GetProperTime()/CLHEP::picosecond;s.mass[i]=t->GetDynamicParticle()->GetMass()/CLHEP::MeV;s.pdg[i]=t->GetDefinition()->GetPDGEncoding();
  observe(i,s.x[i],s.p[i],s.mass[i],s.pdg[i]);++s.steps;
  if(s.snapshots==0 || s.steps%2000==0)snapshot("transport");
}
inline void vertex(const G4Track& parent,const V& daughter,const std::vector<V>& emitted){
  if(!enabled())return;
  event();auto& s=state();const int i=parent.GetDefinition()->GetPDGEncoding()>0?0:1;
  if(!s.seen[i])throw std::runtime_error("Noether: weak vertex before initial endpoint record");
  move(i,position(parent.GetPosition(),parent.GetGlobalTime()),parent.GetDynamicParticle()->GetCharge()/CLHEP::eplus);
  s.p[i]=momentum(parent.GetDynamicParticle()->Get4Momentum());s.properPs[i]=parent.GetProperTime()/CLHEP::picosecond;
  s.mass[i]=parent.GetDynamicParticle()->GetMass()/CLHEP::MeV;s.pdg[i]=parent.GetDefinition()->GetPDGEncoding();
  observe(i,s.x[i],s.p[i],s.mass[i],s.pdg[i]);snapshot("before_vertex");
  for(const auto& p:emitted){s.flux.p=add(s.flux.p,p);s.flux.j=add(s.flux.j,wedge(s.x[i],p));}
  s.p[i]=daughter;s.skipTrack[i]=parent.GetTrackID();s.pdg[i]=i==0?9900017:-9900017;
  s.mass[i]=static_cast<double>(std::sqrt(daughter[0]*daughter[0]-daughter[1]*daughter[1]-daughter[2]*daughter[2]-daughter[3]*daughter[3]));
  observe(i,s.x[i],s.p[i],s.mass[i],s.pdg[i]);++s.vertices;snapshot("after_vertex");
}
inline void finish(){if(enabled()){
  event();snapshot("end_event");auto& s=state();
  if(s.encounterOutput.is_open()){s.encounterOutput<<"{\"kind\":\"end_event\",\"returns\":"<<s.encounters.returns<<",\"processed_intervals\":"<<s.encounters.processed<<",\"buffer_high_water\":"<<s.encounters.highWater<<"}\n";s.encounterOutput.flush();if(!s.encounterOutput)throw std::runtime_error("encounter output finish failed");}
}}
}
#endif
