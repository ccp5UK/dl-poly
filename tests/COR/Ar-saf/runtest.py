#!/usr/bin/env python3
import dlpoly
import numpy as np

def actual(**kwargs):

    dl = dlpoly.DLPoly(exe=kwargs["exe"],
                   control=kwargs["workdir"]+"Ar.control",
                   config=kwargs["workdir"]+"Ar.config",
                   field=kwargs["workdir"]+"Ar.field",
                   workdir=kwargs["workdir"]+"output")

    dl.run(numProcs=kwargs["np"], mpi=kwargs["mpi"])

    dl.load_correlations()

    ss = np.array(dl.correlations.components[0]['stress_xy-stress_xy'])+np.array(dl.correlations.components[0]['stress_yz-stress_yz'])+np.array(dl.correlations.components[0]['stress_zx-stress_zx'])

    ss = np.array(ss)/ss[0]

    return ss

def expected(**kwargs):
    
    cor = dlpoly.correlations.Correlations(source=kwargs["workdir"]+"Ar.cor")

    ss = np.array(cor.components[0]['stress_xy-stress_xy'])+np.array(cor.components[0]['stress_yz-stress_yz'])+np.array(cor.components[0]['stress_zx-stress_zx'])

    return np.array(ss)/ss[0]
