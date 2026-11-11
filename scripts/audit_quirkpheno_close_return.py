"""Censor unproven first-return intervals from complete QuirkPheno step ROOT files."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from quirksim.close_return import monitor


def audit(path, start_ps=5., arm_radius_mm=.1, lambda_ev=774.):
    import uproot
    if lambda_ev <= 0: raise ValueError('lambda must be positive')
    with uproot.open(path) as root:
        events=root['events'].arrays(['q_steps'],library='ak')
        if len(events)!=1: raise ValueError('only single-event pilot files are supported')
        counts=np.asarray(events['q_steps'][0],int)
        if counts.shape!=(2,): raise ValueError('expected two endpoint step counts')
        a=root['steps'].arrays(['evt','q','t','x','y','z'],library='np')
        if len(np.unique(a['evt']))!=1: raise ValueError('steps must contain exactly one event')
        endpoints=[]
        for q in (0,1):
            mask=a['q']==q
            if np.count_nonzero(mask)!=counts[q]:
                raise ValueError('incomplete steps: require stepPrescale=1')
            endpoints.append((a['t'][mask]*1000.,np.column_stack([a[k][mask] for k in ('x','y','z')])))
        result=monitor(*endpoints[0],*endpoints[1],start_ps=start_ps,
                       arm_radius_mm=arm_radius_mm,match_radius_mm=8*1.973269804e-4/lambda_ev)
    result.update(file=str(path),lambda_ev=lambda_ev,matching_radius_definition='8 hbar*c / Lambda')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root_files',nargs='+',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--start-ps',type=float,default=5.)
    parser.add_argument('--arm-radius-mm',type=float,default=.1)
    parser.add_argument('--lambda-ev',type=float,default=774.)
    args=parser.parse_args()
    results=[]
    for path in args.root_files:
        try: results.append(audit(path,args.start_ps,args.arm_radius_mm,args.lambda_ev))
        except (ValueError, KeyError, OSError) as exc:
            results.append({'file':str(path),'status':'INVALID_INPUT','error':str(exc)})
    args.output.write_text(json.dumps({'scope':'CONDITIONAL_POST_RUN_TRANSPORT_CENSOR',
                                      'cases':results},indent=2)+'\n')
    print(json.dumps({'output':str(args.output),'cases':len(results),
                      'statuses':{s:sum(r['status']==s for r in results) for s in sorted({r['status'] for r in results})}},indent=2))
    if any(r['status']=='INVALID_INPUT' for r in results): raise SystemExit(1)


if __name__=='__main__': main()
