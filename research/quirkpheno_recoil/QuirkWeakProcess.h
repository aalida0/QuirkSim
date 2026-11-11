#pragma once
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include <cstdlib>
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
}
class QuirkWeakProcess : public G4VDiscreteProcess {
public:
  explicit QuirkWeakProcess(const G4ParticleDefinition* neutral);
  G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition*) override;
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;
protected:
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) override { return DBL_MAX; }
private:
  const G4ParticleDefinition* m_neutral;
  G4ParticleChange m_change;
};
