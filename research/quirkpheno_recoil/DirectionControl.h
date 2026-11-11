#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include "TwoBodyKinematics.h"

namespace QPWeak {
inline std::array<double,3> controlledDirection(const std::string& mode, double sign,
    std::array<double,3> momentum, std::array<double,3> planned,
    double parentMass=0.,double productionBoostZ=0.) {
  if (mode!="planned" && mode!="collinear" && mode!="collinear_cm") throw std::domain_error("unknown pion direction mode");
  if (sign!=1. && sign!=-1.) throw std::domain_error("collinear sign must be +/-1");
  double p2=0.,n2=0.;
  for(int i=0;i<3;++i) {
    if(!std::isfinite(momentum[i]) || !std::isfinite(planned[i])) throw std::domain_error("nonfinite direction input");
    p2+=momentum[i]*momentum[i];n2+=planned[i]*planned[i];
  }
  if (!(n2>0) || !std::isfinite(n2) || !std::isfinite(p2)) throw std::domain_error("invalid direction norm");
  if(mode=="collinear_cm") {
    if(!(std::isfinite(parentMass)&&parentMass>0&&std::isfinite(productionBoostZ)&&std::abs(productionBoostZ)<1))
      throw std::domain_error("invalid production-CM control boost or mass");
    const double gamma=1/std::sqrt(1-productionBoostZ*productionBoostZ);
    const Four lab{std::sqrt(parentMass*parentMass+p2),momentum[0],momentum[1],momentum[2]};
    const Four toCM{gamma,0,0,-gamma*productionBoostZ},toLab{gamma,0,0,gamma*productionBoostZ};
    const auto cm=boost(lab,toCM,1.);
    const double pcm=std::sqrt(cm[1]*cm[1]+cm[2]*cm[2]+cm[3]*cm[3]);
    Four axis{1,0,0,0};
    for(int i=0;i<3;++i)axis[i+1]=sign*(pcm>1e-14*parentMass ? cm[i+1]/pcm : planned[i]/std::sqrt(n2));
    const auto inLab=boost(boost(axis,cm,parentMass),toLab,1.);
    const Four inverseParent{lab[0],-lab[1],-lab[2],-lab[3]};
    const auto rest=boost(inLab,inverseParent,parentMass);
    const double norm=std::sqrt(rest[1]*rest[1]+rest[2]*rest[2]+rest[3]*rest[3]);
    return {rest[1]/norm,rest[2]/norm,rest[3]/norm};
  }
  auto out=(mode=="collinear" && p2>0) ? momentum : planned;
  const double norm=std::sqrt((mode=="collinear" && p2>0) ? p2 : n2);
  for(auto& v:out) v=v/norm*(mode=="collinear" ? sign : 1.);
  return out;
}
}
