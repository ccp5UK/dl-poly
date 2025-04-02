#!/usr/bin/env python3
import dlpoly
import os
from pathlib import Path

"""
actual and expected functions for testing a STATIS
file generated from a DL_POLY test case. See also
utils/run-test.py
"""

TEST_STATS = ['Total Extended System Energy', 
              'System Temperature', 
              'Configurational Energy', 
              'Short Range Potential Energy', 
              'Electrostatic Energy', 
              'Chemical Bond Energy', 
              'Valence Angle And 3-Body Potential Energy', 
              'Dihedral, Inversion, And 4-Body Potential Energy', 
              'Tethering Energy', 
              'Enthalpy (Total Energy + Pv)', 
              'Rotational Temperature', 
              'Total Virial', 
              'Short-Range Virial', 
              'Electrostatic Virial', 
              'Bond Virial', 
              'Valence Angle And 3-Body Virial', 
              'Constraint Bond Virial', 
              'Tethering Virial', 
              'Volume',
              'Core-Shell Temperature', 
              'Core-Shell Potential Energy', 
              'Core-Shell Virial', 
              'Md Cell Angle Α', 
              'Md Cell Angle Β', 
              'Md Cell Angle Γ', 
              'Pmf Constraint Virial', 
              'Pressure', 
              'External Degree Of Freedom', 
              'stress xx', 
              'stress xy', 
              'stress xz', 
              'stress yx', 
              'stress yy', 
              'stress yz', 
              'stress zx', 
              'stress zy', 
              'stress zz']

def parse_statis(statis, tsteps=(-1, )):
    return  {f"step {tstep if tstep >= 0 else len(statis[0]) + tstep}: {label}": statis[label][tstep]
             for label in TEST_STATS
             for tstep in tsteps}

def actual(exe: Path, workdir: Path, *, numProcs: int = 1, mpi: str = "mpirun -n"):

    dl = dlpoly.DLPoly(exe=exe,
                   control=workdir/"CONTROL",
                   config=workdir/"CONFIG",
                   field=workdir/"FIELD",
                   output=workdir/"OUTPUT",
                   workdir=workdir,
                   vdw_file=workdir/"TABLE",
                   eam_file=workdir/"TABEAM",
                   ang_file=workdir/"TABANG")

    code = dl.run(numProcs=numProcs, mpi=mpi, debug=True)
    if code != 0:
        raise RuntimeError(f"Error code {code} while running DL_POLY")
    dl.load_statis()

    statis = dl.statis
    return parse_statis(statis)

def expected(reference_dir: Path):
    
    statis = dlpoly.statis.Statis(source=reference_dir/"STATIS")

    return parse_statis(statis)
