import hashlib,io,json,subprocess,sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from quirksim.weak_plan import make_plan,write_plan
from af_validate_inclusive_weak import prepare
from audit_inclusive_weak import parse_vertices,local_minima,root_window


def test_native_all_channels_mass_shells_boosts_and_controls(tmp_path):
    include=Path(__file__).resolve().parents[1]/'research/quirkpheno_recoil';src=tmp_path/'check.cc';exe=tmp_path/'check'
    src.write_text('''#include "GeneralWeakPlan.h"
#include "DirectionControl.h"
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
 try {
 if(argc!=2)return 2;
 const auto pair=QPWeak::loadWeakPlan(argv[1]);
 for(const auto& e:pair.endpoint) {
  QPWeak::Four P{std::sqrt(e.M*e.M+1e10),60000.,0.,80000.};
  for(int control=0;control<3;++control) {
   std::array<double,3> axis{};for(size_t j=1;j<e.products.size();++j)for(int k=0;k<3;++k)axis[k]+=e.products[j].p4[k+1];
   const auto target=QPWeak::controlledDirection("collinear_cm",control==2?-1.:1.,{P[1],P[2],P[3]},axis,e.M,.25);
   const auto products=QPWeak::labProducts(e,P,control?&target:nullptr);QPWeak::Four sum{};
   for(const auto& p:products)for(int k=0;k<4;++k)sum[k]+=p.p4[k];
   for(int k=0;k<4;++k)if(std::abs(sum[k]-P[k])>1e-6)return 3;
   for(size_t j=0;j<products.size();++j) {
     const auto& v=products[j].p4;const double p2=v[1]*v[1]+v[2]*v[2]+v[3]*v[3];
     if(std::abs(v[0]-std::sqrt(products[j].mass*products[j].mass+p2))>1e-6)return 4;
   }
  }
 }
 const auto near=QPWeak::unit({1e-8,0.,-1.});
 const auto r=QPWeak::rotate({0.,0.,1.},{0.,0.,1.},near);
 for(int k=0;k<3;++k)if(std::abs(r[k]-near[k])>1e-13)return 5;
 return 0;
 } catch(const std::exception&) {return 9;}
}''')
    subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic','-I',str(include),str(src),'-o',str(exe)],check=True)
    for mass in (50.,70.,225.,400.):
      for channel in ('electron_neutrino','muon_neutrino')+(() if mass==50 else ('pion',)):
       plan=make_plan(mass,31,channel=channel,proper_ps=(50.,100.));path=tmp_path/'weak.plan';write_plan(path,plan);subprocess.run([str(exe),str(path)],check=True)
    bad=tmp_path/'bad.plan';bad.write_text(path.read_text()+' trailing')
    assert subprocess.run([str(exe),str(bad)],stderr=subprocess.DEVNULL).returncode!=0


def test_preparation_keeps_pairs_and_does_not_call_fixture_an_ensemble(tmp_path):
    state=prepare(tmp_path/'pilot');assert len(state['cases'])==40
    assert 'NOT_POPULATION' in state['scope'];assert state['numerical_validation']=='NOT_RUN'
    assert state['physical_survival_probability'] is None
    for c in state['cases']:
      assert hashlib.sha256(Path(c['macro']).read_bytes()).hexdigest()==c['macro_sha256']
      if not c['weak_enabled']:continue
      siblings=[v for v in state['cases'] if v['channel']==c['channel'] and v['weak_enabled']]
      assert all(v['plan']==c['plan'] for v in siblings)
      assert [p['proper_ps'] for p in c['plan']]==[50.,100.]
      assert hashlib.sha256(Path(c['general_plan_file']).read_bytes()).hexdigest()==c['general_plan_sha256']


def test_malformed_vertex_is_rejected():
    vertex={'kind':'decay','schema':2,'parent_mev':[3,0,0,0],'daughter_mev':[2,0,0,0],'products':[{'p4_mev':[1,0,0,0]}],'proper_ps':50.,'threshold_ps':50.}
    hand={'kind':'handoff','history_preserved':True}
    d,h=parse_vertices(['QP_WEAK_JSON '+json.dumps(vertex),'QP_WEAK_JSON '+json.dumps(hand)]);assert len(d)==len(h)==1
    with pytest.raises(ValueError):parse_vertices(['QP_WEAK_JSON '+json.dumps(vertex)])
    vertex['parent_mev'][0]=3.1
    with pytest.raises(ValueError):parse_vertices(['QP_WEAK_JSON '+json.dumps(vertex),'QP_WEAK_JSON '+json.dumps(hand)])


def test_common_time_local_minimum_against_analytic_straight_lines():
    t=np.linspace(0.,10.,1001);a=np.column_stack((t-5,np.full(len(t),.02),np.zeros(len(t))));b=np.column_stack((5-t,np.zeros(len(t)),np.zeros(len(t))))
    result=local_minima([(t,a),(t,b)]);assert len(result)==1
    assert result[0]['time_cm_ps']==pytest.approx(5.);assert result[0]['separation_mm']==pytest.approx(.02)


def test_streaming_root_window_roundtrip(tmp_path):
    import uproot
    beta=.25;g=1/np.sqrt(1-beta*beta);tc=np.linspace(99.,125.,6001);tl=g*tc
    r=.01*np.sin((tc-102.)*2*np.pi/6.);n=len(tc)
    arrays={'q':np.repeat([0,1],n),'t':np.tile(tl/1000,2),'x':np.r_[r,-r],'y':np.r_[np.full(n,2e-5),np.zeros(n)],'z':np.tile(beta*.299792458*tl,2),'px':np.zeros(2*n),'py':np.zeros(2*n),'pz':np.zeros(2*n),'mass':np.full(2*n,225000.),'pdg':np.r_[np.full(n,9900017),np.full(n,-9900017)],'proper':np.tile(tc/1000,2),'charge':np.zeros(2*n)}
    path=tmp_path/'fixture.root'
    with uproot.recreate(path) as root:root['steps']=arrays
    case={'pair_boost_z':beta,'mass_gev':225.1571702763212,'beta_each_cm':.4,'lambda_ev':4642.}
    dec=[{'time_ps':g*t,'position_mm':[0.,0.,beta*.299792458*g*t]} for t in (60.,100.)]
    a,result=root_window(path,case,dec)
    assert result['first_return']['time_cm_ps']==pytest.approx(102.,abs=1e-4)
    assert result['first_return']['separation_mm']==pytest.approx(2e-5,rel=1e-6)
    assert len(a['q'])<2*n
    with pytest.raises(ValueError,match='budget'):root_window(path,case,dec,maximum_samples=1)
