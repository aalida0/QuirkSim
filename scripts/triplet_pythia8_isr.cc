// QuirkSim standalone Pythia8 LO LHE -> QCD ISR/FSR -> endpoint JSONL.
// No infracolor, weak decays, hadronization, or detector.
#include "Pythia8/Pythia.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace Pythia8;
int main(int argc, char** argv) {
  if(argc<3) { std::cerr<<"Usage: triplet_isr input.lhe output.jsonl [seed] [max_events]\n"; return 2; }
  const int seed=argc>3?std::stoi(argv[3]):12345;
  const int maxEvents=argc>4?std::stoi(argv[4]):0;
  std::ofstream output(argv[2]);
  if(!output || seed<=0 || maxEvents<0) return 2;
  output<<std::setprecision(16);
  Pythia py;
  py.particleData.addParticle(9000111,"q1+","q1-",2,3,0,225.15717);
  py.particleData.addParticle(9000112,"q0","q0~",2,0,0,225.0);
  py.particleData.addParticle(9000113,"q2-","q2+",2,-3,0,225.15717);
  for(int id:{9000111,9000112,9000113}){
    py.particleData.mayDecay(id,false);
    py.particleData.isResonance(id,false);
  }
  py.readString("Beams:frameType = 4");
  py.readString(std::string("Beams:LHEF = ")+argv[1]);
  py.readString("Random:setSeed = on");
  py.readString("Random:seed = "+std::to_string(seed));
  py.readString("PartonLevel:ISR = on");
  py.readString("PartonLevel:FSR = on");
  py.readString("PartonLevel:MPI = off");
  py.readString("HadronLevel:all = off");
  py.readString("SpaceShower:QEDshowerByQ = off");
  py.readString("TimeShower:QEDshowerByQ = off");
  py.readString("TimeShower:QEDshowerByL = off");
  py.readString("TimeShower:QEDshowerByOther = off");
  py.readString("Next:numberShowEvent = 0");
  if(!py.init()) { std::cerr<<"Pythia initialization failed\n"; return 3; }
  int accepted=0, failures=0;
  double ptSum=0;
  while(maxEvents==0||accepted<maxEvents) {
    if(!py.next()) {
      if(py.info.atEndOfFile())break;
      if(++failures>20) { std::cerr<<"Too many failed shower events\n";return 4; }
      continue;
    }
    std::vector<int> indices;
    for(int i=0;i<py.event.size();++i) {
      int id=std::abs(py.event[i].id());
      if(py.event[i].isFinal()&&(id==9000111||id==9000112||id==9000113)) indices.push_back(i);
    }
    if(indices.size()!=2) { std::cerr<<"Not exactly two final quirks\n";return 5; }
    auto &a=py.event[indices[0]],&b=py.event[indices[1]];
    if((a.id()>0)==(b.id()>0))return 6;
    double px=a.px()+b.px(),pyy=a.py()+b.py(),pz=a.pz()+b.pz(),en=a.e()+b.e();
    double pairPt=std::hypot(px,pyy);
    ptSum+=pairPt; ++accepted;
    output<<"{\"schema_version\":1,\"stage\":\"PYTHIA8_LO_QCD_ISR_FSR\","
          <<"\"event_index\":"<<accepted<<",\"pair_pt_gev\":"<<pairPt
          <<",\"pair_mass_gev\":"<<std::sqrt(std::max(0.,en*en-px*px-pyy*pyy-pz*pz))
          <<",\"endpoints\":[";
    for(int k=0;k<2;++k){
      auto &p=py.event[indices[k]];
      int id=p.id(),charge=(std::abs(id)==9000111?1:std::abs(id)==9000113?-1:0)*(id>0?1:-1);
      if(k)output<<",";
      output<<"{\"pdg_id\":"<<id<<",\"charge_e\":"<<charge
            <<",\"quirk_number\":"<<(id>0?1:-1)
            <<",\"mass_gev\":"<<p.m()
            <<",\"p4_gev\":["<<p.e()<<","<<p.px()<<","<<p.py()<<","<<p.pz()<<"]}";
    }
    output<<"]}\n";
  }
  py.stat();
  std::cout<<"accepted="<<accepted<<" failures="<<failures
           <<" mean_pair_pt_gev="<<(accepted?ptSum/accepted:0.)<<"\n";
  return accepted>0?0:7;
}
