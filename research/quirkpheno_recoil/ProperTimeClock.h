#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// Consistent units throughout. A sample is (lab time, total energy).
namespace QPClock {
struct Sample { double time, energy; };
// Integrate m/E dt on the unconsumed trial trajectory. Samples must be pure
// and monotone in lab time as the curve parameter increases.
template<class Sampler>
double properIncrement(double length, double startTime, double mass,
                       Sampler sample, double timeTolerance) {
  if (!(length>=0 && mass>0 && timeTolerance>0 && std::isfinite(length) && std::isfinite(startTime) && std::isfinite(mass) && std::isfinite(timeTolerance))) throw std::domain_error("invalid clock inputs");
  if (length==0) return 0.;
  const auto end=sample(length);
  const double dt=end.time-startTime;
  if (!(std::isfinite(dt) && dt>=0 && std::isfinite(end.energy) && end.energy>0))
    throw std::domain_error("invalid or nonmonotone trial clock");
  if (dt==0) return 0.;
  const double nodes[3]={.5-std::sqrt(15.)/10.,.5,.5+std::sqrt(15.)/10.};
  const double weights[3]={5./18.,4./9.,5./18.};
  double sum=0.;
  for (int n=0;n<3;++n) {
    const double target=startTime+nodes[n]*dt;
    double lo=0,hi=length,tlo=startTime,thi=end.time;
    Sample value=end; bool resolved=false;
    for (int i=0;i<80;++i) {
      double x=lo+(target-tlo)*(hi-lo)/(thi-tlo);
      if (i%4==3 || !(x>lo && x<hi)) x=lo+(hi-lo)/2;
      value=sample(x);
      if (!(std::isfinite(value.time) && std::isfinite(value.energy) && value.energy>0 && value.time>=tlo && value.time<=thi))
        throw std::domain_error("nonmonotone clock sample");
      if (std::abs(value.time-target)<=timeTolerance) { resolved=true; break; }
      if (value.time<target) { lo=x;tlo=value.time; } else { hi=x;thi=value.time; }
    }
    if (!resolved) throw std::runtime_error("lab-time quadrature node unresolved");
    sum+=weights[n]*mass/value.energy;
  }
  return dt*sum;
}
// Bracket a proper-time event to the stated clock tolerance; never consume history.
template<class Clock>
double eventLength(double limit,double target,Clock clock,double tolerance) {
  if (!(limit>0 && target>0 && tolerance>0 && std::isfinite(limit) && std::isfinite(target) && std::isfinite(tolerance))) throw std::domain_error("invalid clock root inputs");
  double lo=0,hi=limit,flo=0,fhi=clock(limit);
  if (!(std::isfinite(fhi) && fhi>=target)) throw std::domain_error("event outside clock bracket");
  for (int i=0;i<100;++i) {
    if (fhi-target<=tolerance) return hi;
    double x=lo+(target-flo)*(hi-lo)/(fhi-flo);
    if (i%4==3 || !(x>lo && x<hi)) x=lo+(hi-lo)/2;
    if (!(x>lo && x<hi)) {
      if (fhi-target<=4*tolerance) return hi;
      throw std::runtime_error("clock root lost representability");
    }
    const double value=clock(x);
    if (!std::isfinite(value)) throw std::domain_error("nonfinite proper clock");
    if (std::abs(value-target)<=tolerance) return x;
    // The quadrature itself has floating-point error at the bracket edges.
    if (value<flo-tolerance || value>fhi+tolerance) throw std::domain_error("nonmonotone proper clock");
    if (value<target) { lo=x;flo=value; } else { hi=x;fhi=value; }
  }
  throw std::runtime_error("proper clock root unresolved");
}
}
