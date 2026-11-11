#include "DirectionControl.h"
#include <cassert>
#include <iostream>

int main() {
  constexpr double M=225.1571702763212,m=225.,mu=.13957039;
  for(double bz:{0.,.25}) for(auto pcm:{std::array<double,3>{98,0,0},std::array<double,3>{5,3,-6},std::array<double,3>{0,0,0}}) {
    const double gg=1/std::sqrt(1-bz*bz);
    QPWeak::Four boostCM{gg,0,0,-gg*bz},boostLab{gg,0,0,gg*bz};
    QPWeak::Four p0{std::sqrt(M*M+pcm[0]*pcm[0]+pcm[1]*pcm[1]+pcm[2]*pcm[2]),pcm[0],pcm[1],pcm[2]};
    auto lab=QPWeak::boost(p0,boostLab,1.);
    for(double sign:{-1.,1.}) {
      auto direction=QPWeak::controlledDirection("collinear_cm",sign,{lab[1],lab[2],lab[3]},{0,1,0},M,bz);
      auto products=QPWeak::pionDecay(M,m,mu,{lab[1],lab[2],lab[3]},direction);
      for(int j=0;j<4;++j)assert(std::abs(products.parent[j]-products.daughter[j]-products.pion[j])<1e-10);
      auto pionCM=QPWeak::boost(products.pion,boostCM,1.);
      auto rest=QPWeak::boost(pionCM,{p0[0],-p0[1],-p0[2],-p0[3]},M);
      double norm=0;for(int j=0;j<3;++j)norm+=pcm[j]*pcm[j];norm=std::sqrt(norm);
      double k=0;for(int j=1;j<4;++j)k+=rest[j]*rest[j];k=std::sqrt(k);
      for(int j=0;j<3;++j)assert(std::abs(rest[j+1]/k-sign*(norm>0 ? pcm[j]/norm : (j==1 ? 1. : 0.)))<1e-10);
    }
  }
  bool threw=false;try{QPWeak::controlledDirection("invalid",1,{1,0,0},{0,1,0});}catch(const std::domain_error&){threw=true;}
  assert(threw);
  std::cout<<"PASS exact CM-collinear pion controls, boosts, both signs and zero-momentum axis\n";
}
