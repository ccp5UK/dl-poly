#!/usr/bin/env bash
#

module purge
module load foss/2022a Python

python3 -m pip uinstall dlpoly-py ruamel.yaml
python3 -m pip install -U git+https://gitlab.com/drFaustroll/dlpoly-py.git ruamel.yaml 

