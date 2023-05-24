#!/usr/bin/env python3
# Test suite for DL_POLY
#   Tests specified in accompanying YAML file test-suite.yml
#
#   Adding tests:
#       Make a path TEST_SET/TEST_CASE/
#       Create a TEST_SET/TEST_CASE/runtest.py with two functions
#           def actual(**kwargs) - compute actual results
#           def expected(**kwargs) - read reference data
#       Add an entry into test-suite.yml
#            tests:
#            - name: "TEST_SET"
#                title: "SOMETHING DESCRIPTIVE"
#                testcases:
#                - name: "TEST_CASE"
#                    title: "SOMETHING MORE SPECIFIC"
#                    type: "regression" # or absolute, or empty (for expected == actual)

from ruamel.yaml import YAML
import numpy as np

import sys
from pathlib import Path

import argparse
import os

parser = argparse.ArgumentParser(description='Test DL_POLY')

parser.add_argument('--exe', type=str, default=None, help='location of DLPOLY.Z executable')
parser.add_argument('--np', type=int, default=1, help='number of processes')
parser.add_argument('--mpi', type=str, default='mpirun -n', help='mpi executable')
parser.add_argument('--verbose', action='store_true', help='verbose output')
parser.add_argument('--save', action='store_true', help='save results')
parser.add_argument('--filepath', type=str, default='test-suite-result.yml', help='save results path')
parser.add_argument('--dir', type=str, default="./", help="location of test-suite.yml")

args = parser.parse_args()

verbose = args.verbose

cwd = os.getcwd()
os.chdir(args.dir)

yaml_parser = YAML()

with open("test-suite.yml", 'rb') as in_file:
    data = yaml_parser.load(in_file)


test_result_string = """\
Total test cases: 0
Test cases passed:
    count: 0
    cases: []
Test cases failed:
    count: 0
    cases: []
"""

results = yaml_parser.load(test_result_string)

for test in data["tests"]:

    if (verbose):
        print("Test set {}".format(test["title"]))

    for case in test["testcases"]:

        results["Total test cases"] += 1

        if (verbose):
            print("\tTest case {}".format(case["title"]))

        testcase = __import__(test["name"]+"."+case["name"]+"."+"runtest", fromlist=[test["name"], case["name"]])

        result = testcase.actual(workdir=test["name"]+"/"+case["name"]+"/", exe=args.exe, np=args.np, mpi=args.mpi)

        expected = testcase.expected(workdir=test["name"]+"/"+case["name"]+"/")

        passed = False

        test_type = "absolute"

        if "type" in case.keys():
            test_type = case["type"]

        if test_type == "regression":

            tol = data["defaults"]["tolerance"]

            if "tolerance" in case.keys():
                tol = case["tolerance"]

            # np.arrays's or int/float will work the same
            mse = np.mean((result-expected)**2)

            passed = mse < tol

            if (verbose):
                print("\t\t result / mse  {} / {}".format(passed, mse))

        else:

            passed = result == expected

            if (verbose):
                print("\t\t result {}".format(passed))

        if (passed):
            results["Test cases passed"]["cases"].append(case["title"])
        else:
            results["Test cases failed"]["cases"].append(case["title"])


results["Test cases passed"]["count"] = len(results["Test cases passed"]["cases"])
results["Test cases failed"]["count"] = len(results["Test cases failed"]["cases"])

p = results["Test cases passed"]["count"]
f = results["Test cases failed"]["count"]
t = results["Total test cases"]

if (p == t):
    print("Status: PASSED\n")
else:
    print("Status: FAILED\n\tPassed {}/{}".format(p, t))

if (args.save):
    yaml_parser.dump(results, Path(args.filepath))
else:
    yaml_parser.dump(results, sys.stdout)

os.chdir(cwd)
