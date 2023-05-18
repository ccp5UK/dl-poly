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

    vv = np.array(dl.correlations.components[0]['velocity_x-velocity_x'])+np.array(dl.correlations.components[0]['velocity_y-velocity_y'])+np.array(dl.correlations.components[0]['velocity_z-velocity_z'])
    vv = np.array(vv)/vv[0]

    return vv

def expected(**kwargs):
    
    cor = dlpoly.correlations.Correlations(source=kwargs["workdir"]+"Ar.cor")

    vv = np.array(cor.components[0]['velocity_x-velocity_x'])+np.array(cor.components[0]['velocity_y-velocity_y'])+np.array(cor.components[0]['velocity_z-velocity_z'])
    vv = np.array(vv)/vv[0]

    return vv
