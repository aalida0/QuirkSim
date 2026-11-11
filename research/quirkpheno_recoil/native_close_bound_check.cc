#include "CloseReturnBound.h"
#include <cassert>
#include <iostream>
#include <limits>
#include <random>

int main() {
  std::mt19937_64 rng(24);
  std::uniform_real_distribution<double> u(-1., 1.), p(0., 1.);
  for (int k = 0; k < 10000; ++k) {
    QPClose::Point a{u(rng),u(rng),u(rng)}, b{u(rng),u(rng),u(rng)};
    QPClose::Point va{.2*u(rng),.2*u(rng),.2*u(rng)}, vb{.2*u(rng),.2*u(rng),.2*u(rng)};
    const double t0=p(rng), t1=t0+p(rng), tp=2.*p(rng)-.5;
    QPClose::Point ae{}, bp{};
    for (int j=0;j<3;++j) { ae[j]=a[j]+va[j]*(t1-t0); bp[j]=b[j]+vb[j]*tp; }
    const double lower=QPClose::lowerBound(a,t0,ae,t1,bp,tp,1.);
    for (int i=0;i<=100;++i) {
      const double t=t0+(t1-t0)*i/100.;
      QPClose::Point aa{},bb{};
      for (int j=0;j<3;++j) { aa[j]=a[j]+va[j]*(t-t0); bb[j]=b[j]+vb[j]*t; }
      assert(lower<=QPClose::norm(aa,bb)+1e-12);
    }
  }
  assert(QPClose::lowerBound({0,0,0},0,{0,0,0},0,{2,0,0},0,1)==2);
  bool threw=false;
  try { QPClose::lowerBound({0,0,0},0,{2,0,0},1,{0,0,0},0,1); }
  catch (const std::domain_error&) { threw=true; }
  assert(threw);
  threw=false;
  try { QPClose::lowerBound({0,0,0},0,{0,0,0},1,{0,0,0},0,std::numeric_limits<double>::quiet_NaN()); }
  catch (const std::domain_error&) { threw=true; }
  assert(threw);
  std::cout << "PASS asynchronous conservative separation bound\n";
}
