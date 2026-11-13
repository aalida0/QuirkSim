#include "EncounterRecorder.h"
#include <iostream>
#include <vector>
using namespace QPEncounter;
Point point(long double t,int i,long double r,int pdg=9900017){Point p;p.t=t;p.x={(i==0?1.L:-1.L)*r/2,(i==0?.0001L:-.0001L),0};p.pdg=i==0?pdg:-pdg;return p;}
int main(int argc,char**){
  if(argc>1){
    Recorder recorder;std::cout.precision(18);std::cout<<'[';bool first=true;
    auto cb=[&](auto e){if(!first)std::cout<<',';first=false;std::cout<<'['<<e.t<<','<<e.miss<<']';};
    for(int i=0;i<2;++i){const long double dt=i==0?.07L:.11L;
      for(int k=0;k*dt<=20;++k){auto t=k*dt;auto p=point(t,i,0);p.x={(i==0?.5L:-.5L)*std::sin(t+(i==0?0.L:.1L)),(i==0?.03L:-.03L)*std::cos(2*t),0};recorder.push(i,p,cb);}}
    std::cout<<"]\n";return 0;
  }
  // Different asynchronous meshes: a straight crossing has exactly known minimum.
  for(bool reverse:{false,true}){
    Recorder r;std::vector<Encounter> got;auto cb=[&](auto e){got.push_back(e);};
    for(int j=0;j<2;++j){int i=reverse?1-j:j;
      for(auto t:(i==0?std::vector<long double>{0,1,2,3,4}:std::vector<long double>{0,.7L,1.5L,2.5L,4}))r.push(i,point(t,i,t-2),cb);
    }
    if(got.size()!=1 || std::abs(got[0].t-2)>1e-12L || std::abs(got[0].miss-.0002L)>1e-12L)throw std::runtime_error("asynchronous crossing");
  }
  // Duplicate-time parent/daughter records prevent interpolation through a decay.
  {Recorder r;std::vector<Encounter> got;auto cb=[&](auto e){got.push_back(e);};
    for(int i=0;i<2;++i){for(auto t:{0.L,1.L,2.L})r.push(i,point(t,i,t-1,17),cb);
      r.push(i,point(2,i,1,9900017),cb);r.push(i,point(3,i,0),cb);r.push(i,point(4,i,-1),cb);}
    if(got.size()!=1 || std::abs(got[0].t-3)>1e-12L)throw std::runtime_error("decay discontinuity");
  }
  {Recorder r(3);auto cb=[](auto){};bool failed=false;
   try{for(int k=0;k<4;++k)r.push(0,point(k,0,k),cb);}catch(const std::runtime_error&){failed=true;}
   if(!failed)throw std::runtime_error("buffer censoring was hidden");}
  {Recorder r;auto cb=[](auto){};r.push(0,point(1,0,0),cb);bool failed=false;
   try{r.push(0,point(0,0,0),cb);}catch(const std::runtime_error&){failed=true;}
   if(!failed)throw std::runtime_error("reversed clock was hidden");}
  std::cout<<"PASS online asynchronous encounters, known minimum, decay jump, buffer and clock failures\n";
}
