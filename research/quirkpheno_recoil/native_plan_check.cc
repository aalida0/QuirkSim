#include "SampledPionPlan.h"
#include <sstream>
#include <cassert>
#include <iostream>
#include <iomanip>
int main(int argc,char** argv) {
  if (argc==2) {
    auto p=loadPionPairPlan(argv[1]);std::cout<<std::setprecision(17)<<'[';
    for(int i=0;i<2;++i) {if(i)std::cout<<',';auto e=p.endpoint[i];std::cout<<'['<<e.properPs<<','<<e.direction[0]<<','<<e.direction[1]<<','<<e.direction[2]<<']';}
    std::cout<<"]\n";return 0;
  }
  std::istringstream valid("QP_PION_PLAN_V1\n0 37.3 0 1 0\n1 76.1 0 0 -1\n");
  auto p=readPionPairPlan(valid);assert(p.endpoint[0].properPs==37.3 && p.endpoint[1].direction[2]==-1);
  for(const char* s : {"BAD 0 1 0 1 0 1 2 0 1 0", "QP_PION_PLAN_V1 0 -1 0 1 0 1 1 0 1 0",
                       "QP_PION_PLAN_V1 1 1 0 1 0 0 1 0 1 0", "QP_PION_PLAN_V1 0 1 0 0 0 1 1 0 1 0",
                       "QP_PION_PLAN_V1 0 1 0 1 0 1 1 0 1 0 extra"}) {
    bool failed=false;try {std::istringstream input(s);readPionPairPlan(input);} catch(...) {failed=true;}assert(failed);
  }
  std::cout<<"PASS sampled-pion plan parsing, endpoint order, lifetimes, unit directions and malformed data\n";
}
