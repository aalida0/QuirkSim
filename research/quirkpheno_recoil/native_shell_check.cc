#include "StableMassShell.h"
#include <array>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <iostream>
using V=std::array<double,4>;
double dot(V a,V b){return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];}
int main(int argc,char**){
 if(argc>1){std::cout<<QPTransport::spatialShellEnabled()<<"\n";return 0;}
 std::array<double,2> maxima{};
 for(int mode=0;mode<2;++mode){double m=225000,T=109000;V p{std::sqrt(m*m+1e10+36e8),100000,0,60000};std::array<long double,4> loss{};
 for(int i=0;i<2000000;++i){V a{1e-6,1e-6*(i%20000<10000?1:-1),0,0},u=p,b{};for(int k=0;k<4;++k)u[k]+=T*a[k];
 double c=2*dot(a,u)/dot(u,u);for(int k=0;k<4;++k)b[k]=c*u[k]-a[k];
 V next=p;for(int k=0;k<4;++k)next[k]+=T*(a[k]-b[k]);
 if(mode==0){double c=m/std::sqrt(dot(next,next));for(auto& v:next)v*=c;}
 else next[0]=QPTransport::shellEnergy(m,next[1],next[2],next[3]);
 // Geant4 stores energy/direction and reconstructs momentum each step.
 double mag=std::hypot(next[1],std::hypot(next[2],next[3])),kin=next[0]-m,dp=std::sqrt(kin*(kin+2*m));
 for(int k=1;k<4;++k)next[k]=dp*(next[k]/mag);
 for(int k=0;k<4;++k)loss[k]+=(long double)next[k]-p[k]-(long double)T*((long double)a[k]-b[k]);
 p=next;
 }
 for(auto x:loss)maxima[mode]=std::max(maxima[mode],std::abs(double(x)));
 std::cout<<mode;for(auto x:loss)std::cout<<' '<<double(x);std::cout<<'\n';
 }
 if(maxima[1]>1e-6 || maxima[0]<100*maxima[1])throw std::runtime_error("spatial impulse precision regression");
}
