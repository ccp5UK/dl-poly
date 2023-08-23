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

    components = dl.correlations.components[0]
    if ('v_x-v_x' in components.keys()):
        vv = np.array(components['v_x-v_x'])+np.array(components['v_y-v_y'])+np.array(components['v_z-v_z'])
    else:
        vv = np.array(components['velocity_x-velocity_x'])+np.array(components['velocity_y-velocity_y'])+np.array(components['velocity_z-velocity_z'])
    
    vv = np.array(vv)/vv[0]

    return vv

def expected(**kwargs):
    
    cor = dlpoly.correlations.Correlations(source=kwargs["workdir"]+"Ar.cor")

    vv = np.array(cor.components[0]['velocity_x-velocity_x'])+np.array(cor.components[0]['velocity_y-velocity_y'])+np.array(cor.components[0]['velocity_z-velocity_z'])
    vv = np.array(vv)/vv[0]

    return vv
