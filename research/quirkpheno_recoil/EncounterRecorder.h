#ifndef QP_ENCOUNTER_RECORDER_H
#define QP_ENCOUNTER_RECORDER_H
// Online common-frame piecewise-linear geometry. No simultaneous-step assumption.
#include <array>
#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <limits>
#include <stdexcept>
namespace QPEncounter {
struct Point {long double t{};std::array<long double,3>x{},p{};long double e{},mass{};int pdg{};};
struct Encounter {long double t{},miss{};std::array<Point,2> endpoint;};
class Recorder {
  std::array<std::deque<Point>,2> q;
  long double done=-std::numeric_limits<long double>::infinity();
  long double lastDerivative=0,lastEnd=-std::numeric_limits<long double>::infinity(),lastReturn=-std::numeric_limits<long double>::infinity();
  bool previousNeutral=false;
  std::size_t limit;
  static Point lerp(const Point& a,const Point& b,long double t){
    if(b.t<=a.t)throw std::runtime_error("zero encounter interpolation interval");
    const auto f=(t-a.t)/(b.t-a.t);Point p=a;p.t=t;
    for(int k=0;k<3;++k){p.x[k]=a.x[k]+f*(b.x[k]-a.x[k]);p.p[k]=a.p[k]+f*(b.p[k]-a.p[k]);}
    p.e=a.e+f*(b.e-a.e);return p;
  }
public:
  std::size_t processed=0,returns=0,highWater=0;
  explicit Recorder(std::size_t maximum=200000):limit(maximum){}
  using Callback=std::function<void(const Encounter&)>;
  void push(int i,const Point& p,const Callback& emit){
    if(i<0 || i>1 || !std::isfinite(p.t))throw std::runtime_error("invalid encounter endpoint/time");
    for(auto v:p.x)if(!std::isfinite(v))throw std::runtime_error("nonfinite encounter position");
    if(!q[i].empty() && p.t<q[i].back().t)throw std::runtime_error("time-reversed encounter stream");
    q[i].push_back(p);highWater=std::max(highWater,q[i].size());
    if(q[i].size()>limit)throw std::runtime_error("asynchronous encounter buffer exceeded; do not downsample");
    while(q[0].size()>1 && q[1].size()>1){
      for(auto& a:q)while(a.size()>1 && a[1].t<=done)a.pop_front();
      if(q[0].size()<2 || q[1].size()<2)break;
      // Duplicate events explicitly delimit a discontinuous decay momentum/PDG.
      bool duplicate=false;
      for(auto& a:q)if(a[0].t==a[1].t){a.pop_front();duplicate=true;}
      if(duplicate)continue;
      const auto lo=std::max({q[0][0].t,q[1][0].t,done});
      const auto hi=std::min(q[0][1].t,q[1][1].t);
      if(hi<=lo){
        bool moved=false;
        for(auto& a:q)if(a[1].t<=lo){a.pop_front();moved=true;}
        if(!moved)break;
        continue;
      }
      std::array<Point,2>a{lerp(q[0][0],q[0][1],lo),lerp(q[1][0],q[1][1],lo)};
      std::array<Point,2>b{lerp(q[0][0],q[0][1],hi),lerp(q[1][0],q[1][1],hi)};
      std::array<long double,3>r{},v{};long double rv=0,vv=0,r2=0;
      for(int k=0;k<3;++k){r[k]=a[0].x[k]-a[1].x[k];v[k]=(b[0].x[k]-b[1].x[k]-r[k])/(hi-lo);rv+=r[k]*v[k];vv+=v[k]*v[k];r2+=r[k]*r[k];}
      const bool neutral=std::abs(q[0][0].pdg)==9900017 && std::abs(q[1][0].pdg)==9900017;
      auto save=[&](long double t){
        if(t-lastReturn<=1e-9L)return;
        Encounter e;e.t=t;long double d=0;
        for(int i=0;i<2;++i)e.endpoint[i]=lerp(q[i][0],q[i][1],t);
        for(int k=0;k<3;++k)d+=std::pow(e.endpoint[0].x[k]-e.endpoint[1].x[k],2);
        e.miss=std::sqrt(d);lastReturn=t;++returns;emit(e);
      };
      if(neutral){
        if(previousNeutral && lastEnd==lo && lastDerivative<=0 && rv>=0 && (lastDerivative<0 || rv>0))save(lo);
        if(vv>0 && rv<0 && rv+vv*(hi-lo)>=0)save(lo-rv/vv);
      }
      lastDerivative=rv+vv*(hi-lo);lastEnd=hi;previousNeutral=neutral;done=hi;++processed;
    }
  }
};
}
#endif
