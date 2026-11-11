#pragma once
#include <array>
#include <cmath>
#include <fstream>
#include <istream>
#include <stdexcept>
#include <string>
struct PionEndpointPlan { double properPs; std::array<double,3> direction; };
struct PionPairPlan { std::array<PionEndpointPlan,2> endpoint; };
inline PionPairPlan readPionPairPlan(std::istream& stream) {
  std::string magic;
  if (!(stream>>magic) || magic!="QP_PION_PLAN_V1") throw std::runtime_error("invalid pion-plan header");
  PionPairPlan plan{};
  for(int i=0;i<2;++i) {
    int slot=-1;auto& p=plan.endpoint[i];
    if (!(stream>>slot>>p.properPs>>p.direction[0]>>p.direction[1]>>p.direction[2]) || slot!=i)
      throw std::runtime_error("missing or reordered pion-plan endpoint");
    double norm=0;for(double x:p.direction) { if(!std::isfinite(x)) throw std::runtime_error("nonfinite direction"); norm+=x*x; }
    if (!(std::isfinite(p.properPs) && p.properPs>0 && std::abs(norm-1)<1e-10))
      throw std::runtime_error("invalid lifetime or nonunit pion direction");
  }
  std::string extra;if(stream>>extra) throw std::runtime_error("trailing pion-plan data");
  return plan;
}
inline PionPairPlan loadPionPairPlan(const std::string& path) {
  std::ifstream stream(path);if(!stream) throw std::runtime_error("cannot open sampled pion plan");
  return readPionPairPlan(stream);
}
