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

    ss = np.array(dl.correlations.components[0]['heat_flux_x-heat_flux_x'])+np.array(dl.correlations.components[0]['heat_flux_y-heat_flux_y'])+np.array(dl.correlations.components[0]['heat_flux_z-heat_flux_z'])

    ss = np.array(ss)/ss[0]

    tc = np.nan

    for d in dl.correlations.derived:
        if 'thermal-conductivity' in d.keys():
            tc = d['thermal-conductivity']['value']

    return {"correlation" : ss, "thermal-conductivity" : tc}

def expected(**kwargs):
    
    cor = dlpoly.correlations.Correlations(source=kwargs["workdir"]+"Ar.cor")

    ss = np.array(cor.components[0]['heat_flux_x-heat_flux_x'])+np.array(cor.components[0]['heat_flux_y-heat_flux_y'])+np.array(cor.components[0]['heat_flux_z-heat_flux_z'])

    for d in cor.derived:
        if 'thermal-conductivity' in d.keys():
            tc = d['thermal-conductivity']['value']

    return {"correlation" : np.array(ss)/ss[0], "thermal-conductivity" : tc}
