#include "NoetherCharges.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace QPNoether;
long double dot(const V& a,const V& b){return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];}
V reflect(const V& a,const V& axis){return add(scale(axis,2*dot(a,axis)/dot(axis,axis)),scale(a,-1));}
V ray(std::mt19937_64& g){std::uniform_real_distribution<long double> u(-1,1);V a{0,u(g),u(g),u(g)};a=scale(a,.001L/std::sqrt(-dot(a,a)));a[0]=.001L;return a;}
long double error(const Charge& a,const Charge& b){long double e=0;for(int k=0;k<4;++k)e=std::max(e,std::abs(a.p[k]-b.p[k]));for(int k=0;k<16;++k)e=std::max(e,std::abs(a.j[k]-b.j[k]));return e;}
int main(){
  // Independent straight, static-string integral: E=T L, M^{01}=-T L^2/2.
  {std::array<V,2> x{V{0,0,0,0},V{0,2,0,0}},p{};
   std::array<std::vector<V>,2> rays{std::vector<V>{V{1,1,0,0}},std::vector<V>{V{1,-1,0,0}}};
   auto q=total(x,p,rays,3);
   if(q.p[0]!=6 || q.p[1]!=0 || q.j[1]!=-6 || q.j[4]!=6)throw std::runtime_error("static string stress integral");}
  std::mt19937_64 rng(83117);long double maximum=0;
  for(int trial=0;trial<200;++trial){
    const long double tension=109000,m=225000;
    std::array<std::vector<V>,2> rays;
    for(auto& q:rays)for(int k=0;k<5;++k)q.push_back(ray(rng));
    std::array<V,2> x{V{20,3,-7,2},V{20,3,-7,2}};
    for(auto a:rays[0])x[1]=add(x[1],a);
    for(auto a:rays[1])x[1]=add(x[1],scale(a,-1));
    std::array<V,2> p{V{std::sqrt(m*m+100000.L*100000),100000,0,0},V{std::sqrt(m*m+80000.L*80000),0,-80000,0}};
    const auto initial=total(x,p,rays,tension);
    // Alternate endpoint steps; consume partial and full queue vectors.
    for(int step=0;step<100;++step){int i=step%2;auto a=rays[i].front();
      if(step%3==0){a=scale(a,.37L);rays[i][0]=scale(rays[i][0],.63L);}else rays[i].erase(rays[i].begin());
      auto b=reflect(a,add(p[i],scale(a,tension)));
      p[i]=add(p[i],scale(add(a,scale(b,-1)),tension));x[i]=add(x[i],add(a,b));rays[1-i].push_back(b);
      auto now=total(x,p,rays,tension);maximum=std::max(maximum,error(initial,now));
      for(auto z:closure(x,rays))if(std::abs(z)>1e-12L)throw std::runtime_error("contour did not close");
      if(std::abs(dot(p[i],p[i])-m*m)>1e-5L)throw std::runtime_error("endpoint mass shell");
    }
    // A pointlike decay replaces one endpoint by daughter + emitted products.
    auto before=total(x,p,rays,tension);V emitted{140,30,40,0};p[0]=add(p[0],scale(emitted,-1));
    auto after=total(x,p,rays,tension);after.p=add(after.p,emitted);after.j=add(after.j,wedge(x[0],emitted));
    if(error(before,after)>1e-8L)throw std::runtime_error("vertex flux");
    // A missing queue must be detected; endpoint-only is not a valid shortcut.
    auto missing=rays;missing[0].clear();if(error(total(x,p,rays,tension),total(x,p,missing,tension))<1)throw std::runtime_error("missing string not detected");
  }
  if(maximum>1e-7L)throw std::runtime_error("causal charges drifted");
  std::cout<<"PASS 200 causal histories, 20000 reflections, partial queues, vertex flux; max native drift "<<double(maximum)<<'\n';
}
