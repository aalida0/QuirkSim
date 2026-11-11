#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

// Conservative separation bound during one endpoint step, using one
// asynchronous partner observation. Positions and c*time use the same units.
// Requires both endpoint worldlines to be continuous and speed <= c.
// This is a mathematical bound in exact arithmetic, not a floating-point
// certificate. The helper is not wired to the Geant4 track lifecycle.
namespace QPClose {
using Point=std::array<double,3>;
inline double norm(const Point& a,const Point& b) {
  double s=0;for(int i=0;i<3;++i)s+=(a[i]-b[i])*(a[i]-b[i]);return std::sqrt(s);
}
inline double lowerBound(Point start,double t0,Point end,double t1,
                         Point partner,double tp,double c) {
  if (!(std::isfinite(t0)&&std::isfinite(t1)&&std::isfinite(tp)&&std::isfinite(c)&&c>0&&t1>=t0))
    throw std::domain_error("invalid close-return times or speed");
  for(auto p:{start,end,partner}) for(double x:p) if(!std::isfinite(x))throw std::domain_error("nonfinite position");
  const double dt=t1-t0;
  if(norm(start,end)>c*dt+1e-12*std::max(1.,c*dt))throw std::domain_error("superluminal endpoint step");
  const double partnerTravel=c*std::max(std::abs(t0-tp),std::abs(t1-tp));
  return std::max({0.,norm(start,partner)-c*dt-partnerTravel,
                      norm(end,partner)-c*dt-partnerTravel});
}
}
