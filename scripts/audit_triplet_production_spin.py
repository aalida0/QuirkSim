"""Specified DY source density and local pion spin maps from recorded vertices."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import tarfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from quirksim.production_spin import dy_vector_density,singlet_fraction,pion_spin_map,canonical_boost_rotation,rotate_pair
from quirksim.physics import benchmark


def calculate(archive):
    mass=benchmark('triplet',225.)['charged_mass_gev']
    p=np.array([mass*.4/math.sqrt(1-.4**2),0,0])
    rho=dy_vector_density(p,mass)
    maps=[canonical_boost_rotation(pp,mass,[0,0,.25]) for pp in (p,-p)]
    labrho=rotate_pair(rho,*maps)
    result={'scope':'SPECIFIED_VECTOR_CURRENT_DY_AND_LOCAL_DECAY_AMPLITUDES',
            'production_model':{'mass_gev':mass,'pair_cm_beta_each':.4,'pair_lab_boost_z':.25,
                                'incoming_quarks':'massless; transverse current; equal helicity weights',
                                'final_heavy_current':'pure vector; Dirac Y=0 triplet charged pair',
                                'included_in_AF_gun_generator':False},
            'production_pair_cm_density_real':rho.real.tolist(),
            'production_pair_cm_density_imag':rho.imag.tolist(),
            'production_pair_cm_singlet_fraction':singlet_fraction(rho),
            'lab_canonical_singlet_fraction':singlet_fraction(labrho),
            'frame_warning':'Lab spin-only singlet fraction is not the pair-CM annihilation occupation',
            'source_archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
            'local_decay_maps':[],'spin_density_at_first_return':None,
            'physical_annihilation_probability':None,
            'missing':['continuous spin transport under string and EM forces',
                       'spin-dependent radiation and capture','finite-energy spin-one annihilation',
                       'production wavepacket and spatial matching']}
    with tarfile.open(archive) as tar:
        manifest=json.load(tar.extractfile('manifest.json'))
        for case in manifest['cases']:
            log=tar.extractfile(case['name']+'.log').read().decode()
            for line in log.splitlines():
                if 'QP_WEAK_JSON ' not in line:continue
                event=json.loads(line.split('QP_WEAK_JSON ',1)[1])
                if event['kind']!='decay':continue
                parent=np.array(event['parent_mev'])/1000;daughter=np.array(event['daughter_mev'])/1000
                rotation=pion_spin_map(parent[1:],event['parent_mass_mev']/1000,
                                        daughter[1:],event['daughter_mass_mev']/1000)
                angle=2*math.atan2(float(np.linalg.norm([rotation[0,1],rotation[0,0].imag])),
                                   float(np.trace(rotation).real/2))
                result['local_decay_maps'].append({'case':case['name'],'time_ps':event['time_ps'],
                    'parent_pdg':event['parent_pdg'],'canonical_rotation_angle_rad':angle,
                    'rotation_real':rotation.real.tolist(),'rotation_imag':rotation.imag.tolist(),
                    'map_scope':'LOCAL_VERTEX_ONLY; NO_CONTINUOUS_SPIN_TRANSPORT'})
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=calculate(a.archive);a.output.write_text(json.dumps(r,indent=2)+'\n')
    print('Production pair-CM singlet fraction:',r['production_pair_cm_singlet_fraction'])
    print('Lab canonical spin-only singlet fraction:',r['lab_canonical_singlet_fraction'])
    print('Recorded pion vertices:',len(r['local_decay_maps']))
    print('Max local canonical rotation:',max(x['canonical_rotation_angle_rad'] for x in r['local_decay_maps']))
