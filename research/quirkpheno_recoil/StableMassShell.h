#ifndef QP_STABLE_MASS_SHELL_H
#define QP_STABLE_MASS_SHELL_H
#include <cmath>
#include <cstdlib>
#include <cstring>
namespace QPTransport {
// Opt-in until the native two-resolution AF preflight establishes the change.
inline bool spatialShellEnabled(){
  static const bool enabled=[](){const auto* s=std::getenv("QP_SPATIAL_MASS_SHELL");return s && std::strcmp(s,"1")==0;}();
  return enabled;
}
// The exact string reflection already preserves the mass shell. Retain its
// spatial impulse; recompute E without repeatedly rescaling all p components.
inline double shellEnergy(double mass,double px,double py,double pz){
  return std::hypot(mass,std::hypot(px,std::hypot(py,pz)));
}
}
#endif
