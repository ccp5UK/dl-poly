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
    if ('s_xy-s_xy' in components.keys()):
        ss = np.array(components['s_xy-s_xy'])+np.array(components['s_yz-s_yz'])+np.array(components['s_zx-s_zx'])
    else:
        ss = np.array(components['stress_xy-stress_xy'])+np.array(components['stress_yz-stress_yz'])+np.array(components['stress_zx-stress_zx'])

    ss = np.array(ss)/ss[0]
    
    visc = np.nan

    for d in dl.correlations.derived:
        if 'viscosity' in d.keys():
            visc = d['viscosity']['value']

    return {"correlation" : ss, "viscosity" : visc}

def expected(**kwargs):
    
    cor = dlpoly.correlations.Correlations(source=kwargs["workdir"]+"Ar.cor")

    ss = np.array(cor.components[0]['stress_xy-stress_xy'])+np.array(cor.components[0]['stress_yz-stress_yz'])+np.array(cor.components[0]['stress_zx-stress_zx'])

    for d in cor.derived:
        if 'viscosity' in d.keys():
            visc = d['viscosity']['value']

    return {"correlation" : np.array(ss)/ss[0], "viscosity" : visc}
