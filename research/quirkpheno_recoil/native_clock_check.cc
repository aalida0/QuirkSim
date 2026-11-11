#include "ProperTimeClock.h"
#include <cassert>
#include <iostream>
#include <cmath>
int main() {
  auto constant=[](double s) { return QPClock::Sample{10.+s+.1*s*s,5.}; };
  auto clock=[&](double s) { return QPClock::properIncrement(s,10.,3.,constant,1e-14); };
  assert(std::abs(clock(.4)-.6*(.4+.1*.4*.4))<1e-14);
  const double length=QPClock::eventLength(.4,.2,clock,1e-12);
  const double expected=(-1.+std::sqrt(1.+.4*(.2/.6)))/.2;
  assert(std::abs(length-expected)<1e-10 && std::abs(clock(length)-.2)<1e-12);
  // Independent analytic proper time for a relativistic constant-force orbit.
  const double mass=2.,p0=3.,force=.8;
  for (double start : {0.,1.,3.,3.75,4.}) {
    auto sample=[&](double s) { double t=start+s+.1*s*s; return QPClock::Sample{t,std::hypot(mass,p0-force*t)}; };
    for (double limit : {.0001,.001,.01,.02}) {
      const double t1=sample(limit).time;
      const double exact=mass/force*(std::asinh((p0-force*start)/mass)-std::asinh((p0-force*t1)/mass));
      const double value=QPClock::properIncrement(limit,start,mass,sample,1e-14);
      assert(std::abs(value-exact)<5e-14);
      auto local=[&](double s) { return QPClock::properIncrement(s,start,mass,sample,1e-14); };
      const double target=exact*.67;
      const double at=QPClock::eventLength(limit,target,local,1e-12);
      assert(std::abs(local(at)-target)<1.01e-12 && at<=limit);
    }
  }
  bool failed=false;
  try { QPClock::eventLength(.4,1.,clock,1e-12); } catch (...) { failed=true; } assert(failed);
  failed=false;
  try { QPClock::properIncrement(.1,10.,3.,[](double s) { return QPClock::Sample{10.-s,5.}; },1e-14); }
  catch (...) { failed=true; } assert(failed);
  failed=false;
  try { QPClock::properIncrement(.1,0.,3.,[](double s) { return QPClock::Sample{s,-1.}; },1e-14); }
  catch (...) { failed=true; } assert(failed);
  failed=false;
  try { QPClock::properIncrement(.1,10.,std::numeric_limits<double>::infinity(),constant,1e-14); }
  catch (...) { failed=true; } assert(failed);
  failed=false;
  try { QPClock::properIncrement(.1,0.,3.,[](double) { return QPClock::Sample{std::numeric_limits<double>::quiet_NaN(),5.}; },1e-14); }
  catch (...) { failed=true; } assert(failed);
  assert(QPClock::properIncrement(0.,10.,3.,constant,1e-14)==0.);
  std::cout << "PASS native proper-clock quadrature, analytic force reference, event roots and invalid brackets\n";
}
