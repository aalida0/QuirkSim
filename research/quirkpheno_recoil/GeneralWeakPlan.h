#pragma once
#include "TwoBodyKinematics.h"
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
namespace QPWeak {
struct WeakProduct { int pdg{};double mass{};Four p4{}; };
struct WeakEndpoint { double properPs{},M{},m{};std::string channel;std::vector<WeakProduct> products; };
struct WeakPair { std::array<WeakEndpoint,2> endpoint; };
inline WeakPair readWeakPlan(std::istream& in) {
  std::string magic;if(!(in>>magic)||magic!="QP_WEAK_PLAN_V2")throw std::runtime_error("invalid inclusive weak plan header");
  WeakPair plan;
  for(int i=0;i<2;++i) {
    auto& e=plan.endpoint[i];int slot=-1,n=0;
    if(!(in>>slot>>e.properPs>>e.channel>>e.M>>e.m>>n)||slot!=i)throw std::runtime_error("invalid endpoint row");
    if(!(std::isfinite(e.properPs)&&e.properPs>0&&std::isfinite(e.M)&&std::isfinite(e.m)&&e.M>e.m&&e.m>0))throw std::runtime_error("invalid endpoint scales");
    if(e.channel!="pion"&&e.channel!="electron_neutrino"&&e.channel!="muon_neutrino")throw std::runtime_error("unknown weak channel");
    if(n!=(e.channel=="pion"?2:3))throw std::runtime_error("product count mismatch");
    Four sum{};
    for(int j=0;j<n;++j) {
      WeakProduct p;
      if(!(in>>p.pdg>>p.mass>>p.p4[0]>>p.p4[1]>>p.p4[2]>>p.p4[3]))throw std::runtime_error("missing product");
      if(!std::isfinite(p.mass)||p.mass<0)throw std::runtime_error("invalid product mass");
      for(int k=0;k<4;++k) {if(!std::isfinite(p.p4[k]))throw std::runtime_error("nonfinite product");sum[k]+=p.p4[k];}
      if(p.p4[0]<p.mass)throw std::runtime_error("submass energy");
      const double p2=p.p4[1]*p.p4[1]+p.p4[2]*p.p4[2]+p.p4[3]*p.p4[3];
      if(std::abs(p.p4[0]-std::sqrt(p.mass*p.mass+p2))>1e-6)throw std::runtime_error("off-shell product");
      e.products.push_back(p);
    }
    const int sign=i==0?1:-1,lep=e.channel=="electron_neutrino"?11:13;
    if(e.products[0].pdg!=sign*9900017||std::abs(e.products[0].mass-e.m)>1e-6)throw std::runtime_error("neutral mismatch");
    if(e.channel=="pion") {if(e.products[1].pdg!=sign*211)throw std::runtime_error("pion charge mismatch");}
    else if(e.products[1].pdg!=-sign*lep||e.products[2].pdg!=sign*(lep+1))throw std::runtime_error("lepton flavor/charge mismatch");
    for(int k=0;k<4;++k)if(std::abs(sum[k]-(k==0?e.M:0.))>1e-6)throw std::runtime_error("rest conservation failed");
  }
  std::string extra;if(in>>extra)throw std::runtime_error("trailing weak-plan content");return plan;
}
inline WeakPair loadWeakPlan(const std::string& path) {std::ifstream f(path);if(!f)throw std::runtime_error("cannot open inclusive weak plan");return readWeakPlan(f);}
inline std::array<double,3> unit(std::array<double,3> v) {double n=0;for(double x:v)n+=x*x;if(!(n>0&&std::isfinite(n)))throw std::runtime_error("invalid direction");for(auto& x:v)x/=std::sqrt(n);return v;}
inline std::array<double,3> cross(std::array<double,3> a,std::array<double,3> b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline double dot(std::array<double,3> a,std::array<double,3> b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline std::array<double,3> rotate(std::array<double,3> v,std::array<double,3> a,std::array<double,3> b) {
  a=unit(a);b=unit(b);const double c=std::clamp(dot(a,b),-1.,1.);
  auto axis=cross(a,b);const double sine=std::sqrt(dot(axis,axis));
  if(sine<1e-15) {
    if(c>0)return v;
    int k=0;for(int j=1;j<3;++j)if(std::abs(a[j])<std::abs(a[k]))k=j;
    std::array<double,3> e{};e[k]=1.;axis=unit(cross(a,e));
    const double d=dot(axis,v);for(int j=0;j<3;++j)v[j]=2*d*axis[j]-v[j];
    return v;
  }
  for(auto& x:axis)x/=sine;
  const auto crossv=cross(axis,v);const double axial=dot(axis,v);
  for(int j=0;j<3;++j)v[j]=c*v[j]+sine*crossv[j]+(1-c)*axial*axis[j];
  return v;
}
inline std::vector<WeakProduct> labProducts(const WeakEndpoint& e,Four parent,const std::array<double,3>* target=nullptr) {
  std::array<double,3> axis{};for(size_t j=1;j<e.products.size();++j)for(int k=0;k<3;++k)axis[k]+=e.products[j].p4[k+1];
  auto out=e.products;
  for(auto& p:out) {
    if(target) {const auto v=rotate({p.p4[1],p.p4[2],p.p4[3]},axis,*target);for(int k=0;k<3;++k)p.p4[k+1]=v[k];}
    p.p4=boost(p.p4,parent,e.M);
  }
  return out;
}
}
