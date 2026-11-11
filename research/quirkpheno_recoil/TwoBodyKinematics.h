#pragma once
#include <array>
#include <cmath>
#include <stdexcept>

namespace QPWeak {
using Four = std::array<double,4>; // E,px,py,pz; consistent energy units
struct Products { Four parent, daughter, pion; };
inline Four boost(Four v, const Four& parent, double mass) {
  const double b2=(parent[1]*parent[1]+parent[2]*parent[2]+parent[3]*parent[3])/(parent[0]*parent[0]);
  if (b2==0) return v;
  const double gamma=parent[0]/mass;
  const double bp=(parent[1]*v[1]+parent[2]*v[2]+parent[3]*v[3])/parent[0];
  Four out{}; out[0]=gamma*(v[0]+bp);
  for (int i=1;i<4;++i) out[i]=v[i]+((gamma-1)*bp/b2+gamma*v[0])*parent[i]/parent[0];
  return out;
}
inline Products pionDecay(double M, double m, double mu,
                         std::array<double,3> momentum, std::array<double,3> direction) {
  if (!(m>0 && mu>0 && M>m+mu)) throw std::domain_error("closed pion channel");
  double norm=0, p2=0;
  for (int i=0;i<3;++i) { norm+=direction[i]*direction[i]; p2+=momentum[i]*momentum[i]; }
  if (!(norm>0) || !std::isfinite(norm) || !std::isfinite(p2)) throw std::domain_error("invalid decay vector");
  norm=std::sqrt(norm);
  const double k=std::sqrt((M-m-mu)*(M+m+mu)*(M-m+mu)*(M+m-mu))/(2*M);
  Four parent{std::sqrt(M*M+p2),momentum[0],momentum[1],momentum[2]};
  Four d{std::hypot(m,k),0,0,0}, pi{std::hypot(mu,k),0,0,0};
  for (int i=0;i<3;++i) { pi[i+1]=k*direction[i]/norm; d[i+1]=-pi[i+1]; }
  return {parent,boost(d,parent,M),boost(pi,parent,M)};
}
}
