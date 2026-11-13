#ifndef QP_NOETHER_CHARGES_H
#define QP_NOETHER_CHARGES_H
// Classical orbital charges on the causal contour joining the endpoint events.
// X=(ct,x,y,z) [mm], P=(E,px,py,pz) [MeV], tension [MeV/mm].
#include <array>
#include <vector>
namespace QPNoether {
using V=std::array<long double,4>;
using T=std::array<long double,16>;
inline V add(V a,const V& b){for(int k=0;k<4;++k)a[k]+=b[k];return a;}
inline V scale(V a,long double s){for(auto& x:a)x*=s;return a;}
inline T wedge(const V& x,const V& p){T a{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)a[4*i+j]=x[i]*p[j]-x[j]*p[i];return a;}
inline T add(T a,const T& b){for(int k=0;k<16;++k)a[k]+=b[k];return a;}
struct Charge { V p{};T j{}; };
inline Charge queue(const V& endpoint,const std::vector<V>& segments,long double tension){
  Charge q;V x=endpoint;
  for(const auto& s:segments){q.p=add(q.p,scale(s,tension));q.j=add(q.j,wedge(x,scale(s,tension)));x=add(x,s);}
  return q;
}
inline Charge total(const std::array<V,2>& x,const std::array<V,2>& p,
                    const std::array<std::vector<V>,2>& rays,long double tension){
  Charge q;
  for(int i=0;i<2;++i){auto s=queue(x[i],rays[i],tension);q.p=add(q.p,add(p[i],s.p));q.j=add(q.j,add(wedge(x[i],p[i]),s.j));}
  return q;
}
inline V closure(const std::array<V,2>& x,const std::array<std::vector<V>,2>& rays){
  V c=add(x[0],scale(x[1],-1));
  for(const auto& s:rays[0])c=add(c,s);
  for(const auto& s:rays[1])c=add(c,scale(s,-1));
  return c;
}
}
#endif
