#!/usr/bin/env python3
"""Build/run standalone Pythia8 ISR from MG5 LHE, no MG5 Pythia plugin needed."""
import argparse
import gzip
import os
from pathlib import Path
import shlex
import shutil
import subprocess

def flags():
    config=shutil.which("pythia8-config")
    if config:
        return (shlex.split(subprocess.check_output([config,"--cxxflags"],text=True)),
                shlex.split(subprocess.check_output([config,"--libs"],text=True)))
    if shutil.which("pkg-config") and subprocess.run(
            ["pkg-config","--exists","pythia8"],check=False).returncode==0:
        return (shlex.split(subprocess.check_output(["pkg-config","--cflags","pythia8"],text=True)),
                shlex.split(subprocess.check_output(["pkg-config","--libs","pythia8"],text=True)))
    prefix=os.environ.get("PYTHIA8_DIR","")
    if prefix and (Path(prefix)/"include/Pythia8/Pythia.h").exists():
        lib=Path(prefix)/"lib"
        return ([f"-I{prefix}/include"],[f"-L{lib}",f"-Wl,-rpath,{lib}","-lpythia8"])
    raise RuntimeError("Pythia8 headers/libs not found. Set PYTHIA8_DIR to Pythia8 install prefix.")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lhe",type=Path)
    parser.add_argument("--out",type=Path,default=Path("results/triplet_pythia_isr.jsonl"))
    parser.add_argument("--seed",type=int,default=12345)
    parser.add_argument("--nevents",type=int,default=0,help="0 means all LHE events")
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    inp=args.lhe.resolve(strict=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    build=root/"build"
    build.mkdir(exist_ok=True)
    if inp.suffix==".gz":
        flat=build/(inp.stem+"_pythia_input.lhe")
        with gzip.open(inp,"rb") as source,flat.open("wb") as dest:
            shutil.copyfileobj(source,dest)
        inp=flat
    cflags,libs=flags()
    executable=build/"triplet_pythia8_isr"
    cmd=["g++","-std=c++17","-O2","-Wall","-Wextra",*cflags,
         str(root/"scripts"/"triplet_pythia8_isr.cc"),"-o",str(executable),*libs]
    print("Building Pythia8 bridge:",shlex.join(cmd),flush=True)
    subprocess.run(cmd,check=True)
    cmd=[str(executable),str(inp),str(args.out.resolve()),str(args.seed),str(args.nevents)]
    print("Running:",shlex.join(cmd),flush=True)
    subprocess.run(cmd,check=True)
if __name__=="__main__":
    main()
