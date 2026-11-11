#include "EndpointHandoff.h"
#include "TwoBodyKinematics.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
int main() {
  EndpointHandoff h;
  h.begin(7,9900017,11,22);
  assert(h.isParent(7) && !h.isParent(8));
  bool failed=false; try { h.begin(7,9900017,11,22); } catch (...) { failed=true; } assert(failed);
  failed=false; try { h.consume(8,9900017,11,22); } catch (...) { failed=true; } assert(failed && h.pending);
  failed=false; try { h.consume(7,9900017,11,23); } catch (...) { failed=true; } assert(failed && h.pending);
  failed=false; try { h.consume(7,-9900017,11,22); } catch (...) { failed=true; } assert(failed && h.pending);
  h.consume(7,9900017,11,22); assert(!h.pending);
  failed=false; try { h.consume(7,9900017,11,22); } catch (...) { failed=true; } assert(failed);
  h.begin(8,-9900017,22,11); h.consume(8,-9900017,22,11);
  bool first=true; std::cout << std::setprecision(17) << '[';
  for (const auto p : std::vector<std::array<double,3>>{{0,0,0},{98198.1,0,63386.6},{-98198.1,0,63386.6},{100000,-200000,300000}}) {
    for (double sign : {-1.,1.}) {
      const auto k=QPWeak::pionDecay(225157.1702763212,225000.,139.57039,p,{0,sign,0});
      for (int i=0;i<4;++i) assert(std::abs(k.parent[i]-k.daughter[i]-k.pion[i])<1e-6);
      auto m2=[](const QPWeak::Four& v) { return v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3]; };
      assert(std::abs(m2(k.daughter)-225000.*225000.)<.001);
      assert(std::abs(m2(k.pion)-139.57039*139.57039)<1e-5);
      if(!first) std::cout << ',';
      first=false;
      auto vector=[](const QPWeak::Four& v) { std::cout << '[' << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << ']'; };
      std::cout << "{\"parent\":";vector(k.parent);std::cout << ",\"daughter\":";vector(k.daughter);
      std::cout << ",\"pion\":";vector(k.pion);std::cout << ",\"sign\":" << sign << '}';
    }
  }
  failed=false;try { QPWeak::pionDecay(225100.,225000.,139.57039,{0,0,0},{0,1,0}); }catch(...) {failed=true;} assert(failed);
  std::cout << "]\n";
}
