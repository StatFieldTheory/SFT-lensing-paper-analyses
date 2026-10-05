"""Proposed safe wrapper. Historical21node recipe only, not accepted23node build."""
from __future__ import annotations
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path
CANONICAL = Path('/Users/zzhang/projects/angular_statistics/canoes')
MODULE = 'canoes.sachs.sft_input.corr_op.build'
HISTORICAL = ['--omega-m','0.3160919980475834','--h','0.6711','--n-s','0.97','--ell-max','5000','--dense-threshold','500','--n-log-per-decade','16','--lambda-grid','zs5-21pt','--chi-min','50.0','--n-chi-per-pair','256','--accuracy','default','--ell-hi-threshold','800']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def parser():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--canoes-root',type=Path)
 p.add_argument('--python',type=Path,required=True)
 p.add_argument('--pk',type=Path,required=True)
 p.add_argument('--output',type=Path,required=True)
 p.add_argument('--dry-run',action='store_true')
 p.add_argument('--smoke',action='store_true')
 return p

def prepare(args):
 root=(args.canoes_root or Path(os.environ.get('CANOES_ROOT',str(CANONICAL)))).resolve()
 python=args.python.resolve();pk=args.pk.resolve();out=args.output.absolute()
 if out.suffix != '.npz':raise ValueError('Output suffix must be exactly .npz; no implicit np.savez append')
 cli=root/'src/canoes/sachs/sft_input/corr_op/build.py'
 if not root.is_dir() or not cli.is_file():raise ValueError('Invalid canoes root/CLI source')
 if not python.is_file() or not os.access(python,os.X_OK):raise ValueError('Invalid interpreter')
 if not pk.is_file():raise ValueError('Missing PK')
 import numpy as np
 spectrum=np.loadtxt(pk)
 if spectrum.ndim!=2 or spectrum.shape[0]<2 or spectrum.shape[1]!=2 or not np.isfinite(spectrum).all() or np.any(spectrum<=0) or np.any(np.diff(spectrum[:,0])<=0):raise ValueError('PK must be positive finite ordered two-column data')
 # lexists refuses dangling links too; every lifecycle/output path must be fresh.
 targets=[out,out.with_suffix('.meta.json'),Path(str(out)+'.wrapper_log.txt')]+[Path(str(out)+'.wrapper_'+k+'.json') for k in ('started','completed','failed')]
 if any(os.path.lexists(p) for p in targets):raise ValueError('Output/sidecar already exists')
 if not out.parent.is_dir():raise ValueError('Output parent must exist')
 if out.name=='corr_op_table_zs5_21pt_stf_fid_omega03161_h06711.npz':raise ValueError('Active production basename forbidden')
 env=dict(os.environ,PYTHONPATH=str(root/'src'),CANOES_SUPPRESS_METAL_WARNING='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 # Explicit help subprocess is light CLI validation, never a builder invocation.
 help_result=subprocess.run([str(python),'-m',MODULE,'--help'],cwd=root,env=env,capture_output=True,text=True)
 if help_result.returncode:raise ValueError('Interpreter/CLI import or help failed: '+help_result.stderr)
 options=['--pk','--output','--smoke-test'] if args.smoke else ['--pk','--output']+[x for x in HISTORICAL if x.startswith('--')]
 if any(x not in help_result.stdout for x in options):raise ValueError('Required CLI options unavailable')
 cmd=[str(python),'-m',MODULE,'--pk',str(pk)]+(['--smoke-test'] if args.smoke else HISTORICAL)+['--output',str(out)]
 # Pin all imported package source candidates, including background/corr_op dependencies.
 dependencies={str(p):sha(p) for p in sorted((root/'src/canoes').rglob('*.py'))}
 dependencies[str(Path(__file__).resolve())]=sha(__file__);dependencies[str(pk)]=sha(pk)
 versions=subprocess.run([str(python),'-c',"import json,importlib.metadata as m;print(json.dumps({k:m.version(k) for k in ('numpy','scipy','pyccl')}))"],cwd=root,env=env,capture_output=True,text=True)
 if versions.returncode:raise ValueError('Dependency validation failed: '+versions.stderr)
 return dict(command=cmd,cwd=str(root),environment=env,output=str(out),input_source_hashes=dependencies,interpreter_sha256=sha(python),dependency_versions=json.loads(versions.stdout),cli_help_sha256=hashlib.sha256(help_result.stdout.encode()).hexdigest(),recipe='smoke-only' if args.smoke else 'historical zs5-21pt, unchanged settings; not accepted23node candidate',explicit_cli_settings=cmd[3:],defaults='Unspecified defaults belong to pinned CLI/dependency versions; not claimed as historically identical implementation')

def write_new(path,data):
 with Path(path).open('x') as f:json.dump(data,f,indent=2);f.write('\n')
def execute(plan):
 out=Path(plan['output'])
 if out.suffix != '.npz':raise ValueError('Execution output suffix must be exactly .npz')
 started=Path(str(out)+'.wrapper_started.json');record={k:v for k,v in plan.items() if k!='environment'}
 write_new(started,dict(record,started_unix=time.time(),completed=False))
 try:
  if sha(plan['command'][0])!=plan['interpreter_sha256']:raise RuntimeError('Interpreter drift before builder')
  for p,h in plan['input_source_hashes'].items():
   if sha(p)!=h:raise RuntimeError('Input/source drift before builder')
  if os.path.lexists(out) or os.path.lexists(out.with_suffix('.meta.json')):raise RuntimeError('Output became occupied')
  # Reserve destination exclusively. CLI writes only this owned fresh reservation.
  fd=os.open(out,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600);os.close(fd)
  result=subprocess.run(plan['command'],cwd=plan['cwd'],env=plan['environment'],capture_output=True,text=True)
  log=Path(str(out)+'.wrapper_log.txt');write_new(log,dict(stdout=result.stdout,stderr=result.stderr,returncode=result.returncode))
  if result.returncode:raise RuntimeError('Builder returned '+str(result.returncode))
  if not out.is_file() or not out.stat().st_size:raise RuntimeError('Builder output missing/empty')
  if sha(plan['command'][0])!=plan['interpreter_sha256']:raise RuntimeError('Interpreter drift during builder')
  for p,h in plan['input_source_hashes'].items():
   if sha(p)!=h:raise RuntimeError('Input/source drift during builder')
  hashes={str(out):sha(out)};side=out.with_suffix('.meta.json')
  if side.exists():hashes[str(side)]=sha(side)
  write_new(str(out)+'.wrapper_completed.json',dict(record,completed=True,output_hashes=hashes,ended_unix=time.time()))
 except Exception as exc:
  write_new(str(out)+'.wrapper_failed.json',dict(record,completed=False,error=str(exc),ended_unix=time.time(),output_sha256=sha(out) if out.is_file() else None))
  raise

def main():
 plan=prepare(parser().parse_args())
 if '--dry-run' in sys.argv:print(json.dumps({k:v for k,v in plan.items() if k!='environment'},indent=2));return
 execute(plan)
if __name__=='__main__':main()
