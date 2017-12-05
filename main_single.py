# coding: utf-8

###
### Usage:
### mpirun -n 2 python main170616_mpi.py


from sub171111 import *
import multiprocessing
import time
import datetime
import os
import sys
import subprocess

import steps
import steps.model as smod
import steps.utilities.meshio as smeshio
import steps.geom as sgeom
import steps.rng as srng

import steps.mpi
import steps.mpi.solver as psolv
import steps.utilities.geom_decompose as gd
import steps.utilities.metis_support as metis

import math
import random
import logging

import numpy as np
from itertools import count

import steps.solver as ssa_solver

#################################
#################################
def single_simulation(p, mdl, g, NMDAR_VGCC_MolPerArea):
#################################
#################################

    ##
    ## Time stamp
    ##
    ffname = p.DIRECTORY +  "/timestamp.txt"
    with open(ffname,"a") as f:
        d = datetime.datetime.today()
        f.write('Start Time: '+d.strftime('%Y/%m/%d %H:%M:%S')+'\n')

    ## Random seed
    rng = srng.create('mt19937', 16384) 
    rng.initialize(int(time.time()%4294967295)) # The max unsigned long

    ## Simulation Initialization
    sim = ssa_solver.Tetexact(mdl, g.mesh, rng)
    sim.reset()

    ## Initial concentration
    Conc_CaM = 100e-6 # in Molar Unit
    Conc_CB  = 30e-6 # in Molar Unit
    for i in range(g.nCyt_tets):
        sim.setTetConc(g.Cyt_tets[i], 'N0C0', Conc_CaM) # Conc
        sim.setTetConc(g.Cyt_tets[i], 'CB',   Conc_CB) # Conc

    ## Initial concs for surface reactions
    PSDArea   = sim.getPatchArea('Head_PSD') + sim.getPatchArea('Neck_PSD')
    for P_NAME in p.SURF_DOMAIN_NAMES:
        # print P_NAME
        membArea = sim.getPatchArea( P_NAME )
        sim.setPatchAmount(P_NAME,'PA'  , 1.5e-9 * membArea ) # mol/m^2 | umol() 
        sim.setPatchAmount(P_NAME,'Leak', 1.5e-9 * membArea ) # mol/m^2 |


    # 151895.0 molecules .. daitai OK
    #   6000 moleules per 100uM CaM and 0.1fL Spine

    ## Start simulation
    conc   = np.zeros((g.nCyt_tets,p.len_mol))
    tpnts  = np.arange(0.0, p.INT, p.DT)
    ntpnts = tpnts.shape[0]
    print 'Steps:', ntpnts

    ##
    ##
    for j in range(ntpnts):

        ## VGCC / MNDAR stim
        if j == 50:
            Head_PSD_Area = sim.getPatchArea('Head_PSD') 
            Neck_PSD_Area = sim.getPatchArea('Neck_PSD') 
            sim.setPatchAmount('Head_PSD','NR_Glu', NMDAR_VGCC_MolPerArea * Head_PSD_Area)
            sim.setPatchAmount('Neck_PSD','NR_Glu', NMDAR_VGCC_MolPerArea * Neck_PSD_Area)
            tmp1 = sim.getPatchCount('Head_PSD','NR_Glu')
            tmp2 = sim.getPatchCount('Neck_PSD','NR_Glu')
            print 'NMDAR number:   ', tmp1 + tmp2

        ## Simulation
        if steps.mpi.rank == 0:
            d = datetime.datetime.today()
            print 'Step:', j, '; ', str(tpnts[j]*1000), "msec; Simulating",
            print d.strftime('%Y/%m/%d %H:%M:%S')

        ####
        sim.run(tpnts[j])
        ####

        ## Save variables
        if steps.mpi.rank == 0:
            print 'Step:', j, '; ', str(tpnts[j]*1000), "msec; Saving"
        for NAME in p.CYT_MOL_NAME:
            tmpc = []
            for II in g.Cyt_tets :
                tmpc.append( sim.getTetCount(II,NAME) )
            ffname = p.DIRECTORY + "/" + NAME + "_%03i.npy" % j
            np.save(ffname, tmpc)

        ## Save variables for fluxes

        for NAME in p.PSD_DOMAIN_NAMES:
            NMDACa_Head = sim.getPatchSReacExtent(NAME,'NR_O_To_NR_O_Ca')
            ffname = p.DIRECTORY + "/" + NAME + "_NMDACa_%03i.npy" % j
            np.save( ffname, NMDACa_Head )

        for NAME in p.SURF_DOMAIN_NAMES:
            Ca_Remove_Head = sim.getPatchSReacExtent( NAME , 'PA_To_PA_Ca' )
            Ca_Leak_Head   = sim.getPatchSReacExtent( NAME , 'Leak_To_Leak_Ca' )
            ffname = p.DIRECTORY + "/" + NAME + "_CaRemove_%03i.npy" % j
            np.save(ffname, Ca_Remove_Head)
            ffname = p.DIRECTORY + "/" + NAME + "_CaLeak_%03i.npy" % j
            np.save(ffname, Ca_Leak_Head)

        for NAME_MOL in p.CYT_MOL_NAME:
            for NAME_ROI in g.mesh.getAllROINames():
                M_NUM  = sim.getROIAmount(NAME_ROI, NAME_MOL)
                ffname = p.DIRECTORY + "/" + NAME_MOL + "_" + NAME_ROI + "_%03i.npy" % j
                np.save(ffname, M_NUM)

    ##
    ## Time stamp
    ##
    if steps.mpi.rank == 0:
        ffname = p.DIRECTORY +  "/timestamp.txt"
        with open(ffname,"a") as f:
            d = datetime.datetime.today()
            f.write('End Time: '+d.strftime('%Y/%m/%d %H:%M:%S')+'\n')

#################################

if __name__ == '__main__':

    p    = parameters( 0 )
    flux = 8000
    NMDAR_VGCC_MolPerArea = 5.0e-10

    for i in p.IDs:
        p    = parameters( i )
        ID_str = '%s' % i
        # p.DIRECTORY[0] = 'Data1/'+ID_str +'/'
        p.DIRECTORY[0] = 'Data2/'+ID_str +'/'
    
        mdl = gen_model(p, flux)
        g   = gen_geom(p)

        print p.DIRECTORY[0]
        single_simulation(p, mdl, g, NMDAR_VGCC_MolPerArea)


