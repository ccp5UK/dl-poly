#!/usr/bin/env python3

from ruamel.yaml import YAML
import numpy as np

import sys
from pathlib import Path

import argparse
import os

import glob
import subprocess

parser = argparse.ArgumentParser(description='Test DL_POLY')

parser.add_argument('--exe', type=str, default=None, help='location of DLPOLY.Z executable')
parser.add_argument('--np', type=int, default=1, help='number of processes')
parser.add_argument('--mpi', type=str, default='mpirun -n', help='mpi executable')
parser.add_argument('--verbose', action='store_true', help='verbose output')
parser.add_argument('--save', action='store_true', help='save results')
parser.add_argument('--dir', type=str, default="./", help="location of test-suite.yml")

args = parser.parse_args()

verbose = args.verbose

cwd = os.getcwd()
os.chdir(args.dir)

sys.path.append(args.dir)

fname = args.dir.split('/')[-1]

subprocess.run(['tar','-xvJf',f'{args.dir}/{fname}.tar.xz', '--strip-components', '1'])

yaml_parser = YAML()

file = glob.glob(f"{args.dir}/test*.yml")[0]

with open(file, 'rb') as in_file:
    test = yaml_parser.load(in_file)['test']

if (verbose):
    print("Test case {}".format(test["title"]))

test_result_string = """\
Passed: False
Results: {}
"""

results = yaml_parser.load(test_result_string)

testcase = __import__("runtest")

actual = testcase.actual(workdir='./', exe=args.exe, np=args.np, mpi=args.mpi)

expected = testcase.expected(workdir='./')

passed = False

test_type = "absolute"

mse = []

if "type" in test.keys():
    test_type = test["type"]

if (test_type == "absolute"):
    if (a == b):
        passed = True
else:
    tol = test["tolerance"]
    if (isinstance(actual, dict)):
        for i in actual:
            mse.append(np.mean((actual[i]-expected[i])**2))
    else:
        mse.append(np.mean((actual-expected)**2))

    passed = (np.array(mse) < tol).all()

if (isinstance(actual, dict)):
    for (j, i) in enumerate(actual):
        res = {}
        res['Value'] = i
        res['Diff'] = str(mse[j]) if test_type == "regression" else '1'
        results["Results"][i] = res 
else:
    res = {}
    res['Diff'] = str(mse) if test_type == "regression" else '0'
    results["Results"] = res 

results['Passed'] = str(passed)

if (passed):
    print("Status: PASSED\n")
else:
    print("Status: FAILED")

if (args.save):
    yaml_parser.dump(results, Path(args.filepath))
else:
    yaml_parser.dump(results, sys.stdout)

os.chdir(cwd)
