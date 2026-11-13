"""Replayable complete LO weak rest-product plans; MeV in the native file.

Pion: vector-current canonical-rest vertex. Leptons: unpolarized V x V-A.
Forcing a channel/clock is a validation fixture, never an inclusive ensemble.
"""
from pathlib import Path
import numpy as np
from .physics import benchmark,MPI,ME,MMU
from .dynamics import sample_plan
from .decays import sample_decay


def make_plan(mass,seed,channel=None,proper_ps=None):
    b=benchmark('triplet',mass);latent=sample_plan(b['lifetime_ps'],b['widths_gev'],seed);rows=[]
    for i,p in enumerate(latent):
        ch=channel or p['channel'];sign=1 if i==0 else -1
        if b['widths_gev'].get(ch,0)<=0:raise ValueError('closed or unknown forced channel')
        rest=sample_decay(ch,np.zeros(3),b['charged_mass_gev'],mass,np.random.default_rng(p['rest_seed']))
        ids={'daughter':sign*9900017,'pion':sign*211,'charged_lepton':-sign*(11 if ch=='electron_neutrino' else 13),'neutrino':sign*(12 if ch=='electron_neutrino' else 14)}
        masses={'daughter':mass,'pion':MPI,'charged_lepton':ME if ch=='electron_neutrino' else MMU,'neutrino':0.}
        rows.append({'endpoint':i,'proper_ps':float(proper_ps[i] if proper_ps is not None else p['proper_time_ps']),'channel':ch,'M':b['charged_mass_gev'],'m':mass,
            'products':[{'pdg':ids[k],'mass_gev':masses[k],'p4_gev':v.tolist()} for k,v in rest.items()]})
    return rows


def write_plan(path,rows):
    lines=['QP_WEAK_PLAN_V2']
    for i,p in enumerate(rows):
        if p['endpoint']!=i:raise ValueError('ordered pair required')
        lines.append(f"{i} {p['proper_ps']:.17g} {p['channel']} {p['M']*1000:.17g} {p['m']*1000:.17g} {len(p['products'])}")
        for v in p['products']:lines.append(str(v['pdg'])+' '+' '.join(f'{x*1000:.17g}' for x in [v['mass_gev'],*v['p4_gev']]))
    Path(path).write_text('\n'.join(lines)+'\n')
