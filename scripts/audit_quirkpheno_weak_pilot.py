"""Audit deterministic Geant4 handoff diagnostics; no physical survival inference."""
import argparse
import hashlib
import io
import json
import tarfile
import numpy as np
import uproot
from scipy.interpolate import CubicHermiteSpline

def audit(path, plot_path=None):
    results={};curves={};plans={}
    with tarfile.open(path) as archive:
        item=next(m for m in archive.getmembers() if m.name.endswith('/manifest.json'))
        manifest=json.load(archive.extractfile(item));prefix=item.name.rsplit('/',1)[0]+'/'
        for case in manifest['cases']:
            name=case['name'];plans[name]=case
            result={'checks':{},'weak_records':[],'endpoints':[]};results[name]=result;checks=result['checks']
            checks['execution']=case.get('execution',{}).get('returncode')==0 and not case.get('execution',{}).get('timeout',True)
            if prefix+name+'.root' not in archive.getnames():
                checks['root_present']=False;continue
            checks['root_present']=True
            checks['macro_hash']=hashlib.sha256(archive.extractfile(prefix+name+'.mac').read()).hexdigest()==case['macro_sha256']
            log=archive.extractfile(prefix+name+'.log').read().decode(errors='replace')
            checks['no_exception']='*** G4Exception' not in log
            for line in log.splitlines():
                if 'QP_WEAK_JSON ' in line:result['weak_records'].append(json.loads(line.split('QP_WEAK_JSON ',1)[1]))
            decays=[r for r in result['weak_records'] if r['kind']=='decay']
            transfers=[r for r in result['weak_records'] if r['kind']=='handoff']
            expected=int(case['weak_enabled']);checks['decay_count']=len(decays)==expected;checks['handoff_count']=len(transfers)==expected
            for d in decays:
                p=np.array(d['parent_mev']);daughter=np.array(d['daughter_mev']);pion=np.array(d['pion_mev'])
                checks['four_momentum']=bool(max(abs(p-daughter-pion))<1e-6)
                beta=p[1:]/p[0];b2=float(beta@beta);gamma=p[0]/d['parent_mass_mev']
                rest=pion[1:].copy()
                if b2:
                    rest+=((gamma-1)*float(beta@pion[1:])/b2-gamma*pion[0])*beta
                checks['rest_frame_direction']=bool(abs(rest[0])<1e-6 and abs(rest[2])<1e-6 and rest[1]*case['pion_rest_y_sign']>0)
                for label,v in [('parent',p),('daughter',daughter),('pion',pion)]:
                    checks[label+'_on_shell']=bool(abs(v[0]**2-v[1:]@v[1:]-d[label+'_mass_mev']**2)<.01)
                result['observed_parent_mass_gev']=d['parent_mass_mev']/1000
                result['observed_pion_mass_mev']=d['pion_mass_mev']
                result['observed_mass_gap_mev']=d['parent_mass_mev']-d['daughter_mass_mev']
                result['parent_mass_parameter_matches_manifest']=abs(d['parent_mass_mev']-manifest['mass_gev']*1000)<1e-6
                result['proper_threshold_overshoot_ps']=d['proper_ps']-d['threshold_ps']
                checks['threshold']=result['proper_threshold_overshoot_ps']>=0
                checks['endpoint_sign']=(d['daughter_pdg']>0)==(case['endpoint']==0)
                checks['pion_charge']=(d['pion_pdg']>0)==(case['endpoint']==0)
            if expected and len(decays)==1 and len(transfers)==1:
                checks['transfer_identity']=decays[0]['parent_track']==transfers[0]['parent_track'] and decays[0]['daughter_pdg']==transfers[0]['daughter_pdg']
                checks['history_fingerprint']=transfers[0]['history_preserved'] is True
            root=uproot.open(io.BytesIO(archive.extractfile(prefix+name+'.root').read()))
            e=root['events'].arrays(library='ak');a=root['steps'].arrays(library='np')
            checks['one_event']=len(e)==1
            for q in (0,1):
                mask=a['q']==q;tt=a['t'][mask]*1000
                checks[f'end{q}_complete']=int(e['q_reason'][0][q])==6 and float(e['q_te'][0][q])*1000>=manifest['nominal_stop_ps']
                checks[f'end{q}_rows']=len(tt)==int(e['q_steps'][0][q])
                checks[f'end{q}_finite']=bool(np.all(np.isfinite(np.r_[tt,a['px'][mask],a['py'][mask],a['pz'][mask]])))
                checks[f'end{q}_time_order']=bool(np.all(np.diff(tt)>=0))
                neutral=mask & (abs(a['pdg'])==9900017)
                should_neutral=bool(expected and q==case['endpoint'])
                checks[f'end{q}_species']=bool(np.any(neutral))==should_neutral
                if np.any(neutral):
                    checks[f'end{q}_neutral_charge']=bool(np.all(a['charge'][neutral]==0))
                    checks[f'end{q}_neutral_mass']=bool(np.all(abs(a['mass'][neutral]-manifest['neutral_mass_gev']*1000)<1e-6))
                result['endpoints'].append({'q':q,'rows':len(tt),'last_time_ps':float(tt[-1]),'final_pdg':int(a['pdg'][mask][-1])})
                # Fit separate sides of the velocity discontinuity, never a spline across recoil.
                segments=[]
                initial_p=np.array([float(e[k][0][q]) for k in ['q_px0','q_py0','q_pz0']])
                initial_x=np.array([float(e[k][0][q]) for k in ['q_x0','q_y0','q_z0']])
                initial_E=float(e['q_e0'][0][q])
                def curve(indices,start_t,start_x,start_v,end=None):
                    t=np.r_[start_t,a['t'][indices]*1000]
                    x=np.vstack([start_x,np.column_stack([a[k][indices] for k in ['x','y','z']])])
                    pp=np.column_stack([a[k][indices] for k in ['px','py','pz']])
                    v=np.vstack([start_v,.299792458*pp/(a['ekin'][indices]+a['mass'][indices])[:,None]])
                    if end is not None:
                        t=np.r_[t,end[0]];x=np.vstack([x,end[1]]);v=np.vstack([v,end[2]])
                    # Repeated zero-length suspension records contain the same side's state.
                    keep=np.r_[True,np.diff(t)>0];t=t[keep];x=x[keep];v=v[keep]
                    if len(t)>1:return (float(t[0]),float(t[-1]),CubicHermiteSpline(t,x,v,extrapolate=False))
                if should_neutral and len(decays)==1:
                    d=decays[0];td=d['time_ps'];xd=np.array(d['position_mm']);pd=np.array(d['parent_mev']);dd=np.array(d['daughter_mev'])
                    before=mask & (a['handoff']==0) & (abs(a['pdg'])==17) & (a['t']*1000<td)
                    segments.append(curve(before,0,initial_x,.299792458*initial_p/initial_E,(td,xd,.299792458*pd[1:]/pd[0])))
                    segments.append(curve(neutral,td,xd,.299792458*dd[1:]/dd[0]))
                else:
                    segments.append(curve(mask,0,initial_x,.299792458*initial_p/initial_E))
                curves[name,q]=[s for s in segments if s]
            result['local_handoff_checks_pass']=bool(all(checks.values()))
        def value_at(name,q,t):
            for lo,hi,f in curves[name,q]:
                if lo<=t<=hi:return f(t)
            raise ValueError('time outside saved trajectory')
        angle_pairs=[]
        for level in (7,8):
            plus=f'pion_yplus_end0_level{level}';minus=f'pion_yminus_end0_level{level}'
            if not all(results.get(n,{}).get('local_handoff_checks_pass',False) for n in (plus,minus)):continue
            dp=next(r for r in results[plus]['weak_records'] if r['kind']=='decay')
            dm=next(r for r in results[minus]['weak_records'] if r['kind']=='decay')
            same_parent=abs(dp['time_ps']-dm['time_ps'])<1e-10 and np.allclose(dp['parent_mev'],dm['parent_mev'],atol=1e-6,rtol=0)
            stop=manifest['nominal_stop_ps']
            deltas=[(value_at(plus,q,stop)-value_at(minus,q,stop)).tolist() for q in (0,1)]
            td=max(dp['time_ps'],dm['time_ps']);xd=np.array(dp['position_mm'])
            differences=[]
            for t in np.linspace(td,stop,2001):
                xp=value_at(plus,1,t);xm=value_at(minus,1,t)
                # Conservative Euclidean light-cone exclusion, not arrival along the string.
                if .299792458*(t-td)<min(np.linalg.norm(xp-xd),np.linalg.norm(xm-xd))-1e-6:
                    differences.append(float(np.linalg.norm(xp-xm)))
            angle_pairs.append({'level':level,'same_predecay_parent_state_and_time':bool(same_parent),
                                'endpoint_plus_minus_position_delta_at_130ps_mm':deltas,
                                'partner_pre_euclidean_lightcone_max_difference_mm':max(differences) if differences else None,
                                'causal_validation_requires_resolution_and_interpolation_review':True})
        if plot_path is not None:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            grid=np.linspace(0,manifest['nominal_stop_ps'],2601)
            delta={}
            for level in (7,8):
                plus=f'pion_yplus_end0_level{level}';minus=f'pion_yminus_end0_level{level}'
                if not all(results[n]['local_handoff_checks_pass'] for n in (plus,minus)):
                    raise ValueError('Refusing recoil plot for failed handoff checks')
                for q in (0,1):delta[level,q]=np.array([value_at(plus,q,t)[1]-value_at(minus,q,t)[1] for t in grid])
            fig,(ax,bx)=plt.subplots(2,1,figsize=(7.2,6.2),sharex=True,gridspec_kw={'height_ratios':[2,1]},constrained_layout=True)
            for q,color,label in [(0,'#1764ab','decaying endpoint'),(1,'#cf6310','partner')]:
                ax.plot(grid,delta[8,q]*1000,color=color,label=label+' (level 8)')
                ax.plot(grid,delta[7,q]*1000,color=color,ls='--',lw=.9,alpha=.7)
                bx.plot(grid,(delta[8,q]-delta[7,q])*1e6,color=color,label=label)
            decay=next(r for r in results['pion_yplus_end0_level8']['weak_records'] if r['kind']=='decay')
            for axes in (ax,bx):
                axes.axvline(decay['time_ps'],color='.4',ls=':',lw=1)
                axes.grid(alpha=.18);axes.set_xlim(0,130)
            ax.text(decay['time_ps']+2,-.2,'weak decay',color='.35',fontsize=9)
            ax.set_ylabel('Pion +y minus -y position [µm]')
            bx.set_ylabel('Level 8 minus 7 [nm]');bx.set_xlabel('Lab time [ps]')
            ax.legend(frameon=False,fontsize=9);bx.legend(frameon=False,fontsize=8)
            ax.set_title('Deterministic triplet pion recoil, 2 T\nPhysical opposite-direction controls; dashed: level 7',fontsize=11)
            fig.savefig(plot_path)
            plt.close(fig)
        comparison=[]
        for mode in sorted({c['mode'] for c in manifest['cases']}):
            names=[c['name'] for c in sorted(manifest['cases'],key=lambda c:c['level']) if c['mode']==mode]
            if len(names)!=2:continue
            for q in (0,1):
                if not all((n,q) in curves for n in names):continue
                stop=min(manifest['nominal_stop_ps'],*(curves[n,q][-1][1] for n in names))
                grid=np.linspace(0,stop,13001)
                vals=[]
                for n in names:
                    values=np.full((len(grid),3),np.nan)
                    for lo,hi,f in curves[n,q]:
                        take=(grid>=lo)&(grid<=hi);values[take]=f(grid[take])
                    vals.append(values)
                good=np.isfinite(vals[0]).all(axis=1)&np.isfinite(vals[1]).all(axis=1)
                comparison.append({'mode':mode,'endpoint':q,'max_position_difference_mm':float(max(np.linalg.norm(vals[0][good]-vals[1][good],axis=1))),
                                   'sampled_times':int(sum(good)),'stop_ps':float(stop),'interpolation_error_bound_established':False})
    return {'archive_sha256':hashlib.sha256(open(path,'rb').read()).hexdigest(),
            'source_commit':manifest.get('observed_source_commit'),'adapter_sha256':manifest['adapter_sha256'],
            'scope':manifest['scope'],'cases':results,'resolution_comparisons':comparison,
            'paired_direction_comparisons':angle_pairs,
             'local_handoff_checks_pass':all(c.get('local_handoff_checks_pass',False) for c in results.values()),
            'benchmark_mass_parameters_match':all(c.get('parent_mass_parameter_matches_manifest',True) for c in results.values()),
            'causal_propagation_validation':'NOT_ESTABLISHED_BY_LOCAL_HANDOFF_CHECKS',
            'physical_survival_probability':None,
            'note':'Fixed-threshold pion pilot, not sampled lifetimes or a no-recoil control. Stable-versus-decay mixes charge/mass/recoil.'}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive');p.add_argument('--output',required=True);p.add_argument('--plot');a=p.parse_args()
    r=audit(a.archive,a.plot)
    with open(a.output,'w') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
    print('LOCAL_HANDOFF_CHECKS_PASS' if r['local_handoff_checks_pass'] else 'LOCAL_HANDOFF_CHECKS_FAIL')
