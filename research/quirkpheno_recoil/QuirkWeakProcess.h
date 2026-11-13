#pragma once
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include <cstdlib>
#include "quirkcore/SampledPionPlan.h"
#include "quirkcore/GeneralWeakPlan.h"
#include <cfloat>
#include <string>
class G4ParticleDefinition;
namespace QPWeakConfig {
inline bool localizedClock() { const char* p=std::getenv("QP_WEAK_CLOCK_MODE"); return p && std::string(p)=="localized"; }
inline bool enabled() { const char* p=std::getenv("QP_WEAK_ENABLED"); return p && std::string(p)=="1"; }
inline int endpoint() { const char* p=std::getenv("QP_WEAK_ENDPOINT"); return p ? std::stoi(p) : 0; }
inline double tauPs() { const char* p=std::getenv("QP_WEAK_TAU_PS"); return p ? std::stod(p) : 50.; }
inline double neutralGeV() { const char* p=std::getenv("QP_NEUTRAL_MASS_GEV"); return p ? std::stod(p) : 225.; }
inline double directionSign() { const char* p=std::getenv("QP_WEAK_DIR_SIGN"); return p ? std::stod(p) : 1.; }
inline std::string directionMode() { const char* p=std::getenv("QP_WEAK_DIRECTION_MODE"); return p ? p : "planned"; }
inline double collinearSign() { const char* p=std::getenv("QP_WEAK_COLLINEAR_SIGN"); return p ? std::stod(p) : 1.; }
inline double controlBoostZ() { const char* p=std::getenv("QP_WEAK_CONTROL_BOOST_Z"); return p ? std::stod(p) : 0.; }
inline const PionPairPlan* sampledPlan() {
  const char* path=std::getenv("QP_WEAK_PLAN");
  if (!path || !*path) return nullptr;
  static thread_local const PionPairPlan plan=loadPionPairPlan(path);
  return &plan;
}
inline const QPWeak::WeakPair* generalPlan() {
  const char* path=std::getenv("QP_WEAK_GENERAL_PLAN");if(!path||!*path)return nullptr;
  static thread_local const auto plan=QPWeak::loadWeakPlan(path);return &plan;
}
inline double tauPs(int pdg) {
  if(const auto* p=generalPlan())return p->endpoint[pdg>0?0:1].properPs;
  const auto* plan=sampledPlan();return plan ? plan->endpoint[pdg>0 ? 0 : 1].properPs : tauPs();
}
inline std::array<double,3> pionDirection(int pdg) {
  const auto* plan=sampledPlan();return plan ? plan->endpoint[pdg>0 ? 0 : 1].direction : std::array<double,3>{0.,directionSign(),0.};
}
inline bool endpointEnabled(int pdg) {
  return enabled() && std::abs(pdg)==17 && (endpoint()==-1 || (pdg>0 ? 0 : 1)==endpoint());
}
}
class QuirkWeakProcess : public G4VDiscreteProcess {
public:
  explicit QuirkWeakProcess(const G4ParticleDefinition* neutral);
  G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition*) override;
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;
protected:
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) override { return DBL_MAX; }
private:
  G4VParticleChange* generalDecay(const G4Track&);
  const G4ParticleDefinition* m_neutral;
  G4ParticleChange m_change;
};
