##
## python -i main170622_plot.py
##

from sub171111 import *
import math
import sys
import numpy as np
from itertools import count
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import steps
import steps.model as smod
import steps.utilities.meshio as smeshio
import steps.geom as sgeom

import steps.solver as ssolv
import steps.rng as srng
import time


def MolTimeDev(fname, tnum, ID):
    Mol = []
    for j in range(tnum):
        ffname = fname % j
        tmp = sum( np.load(ffname)[ID] )
        Mol.append( tmp.tolist() )
    Mol = np.array( Mol )
    return Mol


def MolTimeFlux(fname, tnum):
    Flux  = []
    prevF = 0
    for j in range(tnum):
        ffname = fname % j
        tmp = np.load(ffname)
        presentF = tmp - prevF
        prevF    = tmp
        Flux.append( presentF.tolist() )
    Flux = np.array( Flux )
    return Flux


def CaFlux( p ):

    ##
    ## Specify target
    ## 
    g   = gen_geom(p)
    tpnts = np.arange(0.0, p.INT, p.DT)
    tpnts = tpnts - tpnts[50]
    tnum  = len( tpnts )

    ##
    ## Load geometry
    ##
    Targ = './Morph/' + p.MorphID
    Head_ID   = np.loadtxt(Targ+'Head_ID.dat') - 1
    Neck_ID   = np.loadtxt(Targ+'Neck_ID.dat') - 1
    Dend_ID   = np.loadtxt(Targ+'Dendrite_ID.dat') - 1
    Dist_Area = np.loadtxt(Targ+'Dist_Area.dat')
    Dist_Area_PSD  = np.loadtxt(Targ+'Face_Loc_Area.dat')

    Head_ID = Head_ID.astype(np.int)
    Neck_ID = Neck_ID.astype(np.int)
    Dend_ID = Dend_ID.astype(np.int)

    ##
    ## Time development of Ca
    ##

    ffname = p.DIRECTORY + "/Ca_%03i.npy"
    Ca_Head = MolTimeDev(ffname, tnum, Head_ID)
    Ca_Neck = MolTimeDev(ffname, tnum, Neck_ID)
    Ca_Dend = MolTimeDev(ffname, tnum, Dend_ID)

    ##
    ## Time development of Ca/CaM
    ##

    CaM_Head = np.zeros( Ca_Head.shape ) # zeros_like(Ca_Head)
    CaM_Neck = np.zeros( Ca_Neck.shape )
    CaM_Dend = np.zeros( Ca_Dend.shape )
    CaM = ['N0C1','N0C2','N1C0','N1C1','N1C2','N2C0','N2C1','N2C2']
    for T_CaM in CaM:
        ffname = p.DIRECTORY + "/"+ T_CaM + "_%03i.npy"
        CaM_Head = CaM_Head + MolTimeDev(ffname, tnum, Head_ID)
        CaM_Neck = CaM_Neck + MolTimeDev(ffname, tnum, Neck_ID)
        CaM_Dend = CaM_Dend + MolTimeDev(ffname, tnum, Dend_ID)

    ###
    ### Ca influx via NMDARs 
    ### Ca uptake-leak through Head/Neck/Dend membranes
    ###

    ffname = p.DIRECTORY + "/Head_PSD_NMDACa_%03i.npy"
    NMDA_Head_Ca = MolTimeFlux(ffname, tnum)
    ffname = p.DIRECTORY + "/Neck_PSD_NMDACa_%03i.npy"
    NMDA_Neck_Ca = MolTimeFlux(ffname, tnum)

    for NAME in p.SURF_DOMAIN_NAMES:
        print 'Removal_'+NAME 
        ffname1 = p.DIRECTORY + "/" + NAME + "_CaRemove_%03i.npy"
        ffname2 = p.DIRECTORY + "/" + NAME + "_CaLeak_%03i.npy"
        ns = globals()
        ns['Removal_'+NAME] = MolTimeFlux(ffname1, tnum) + MolTimeFlux(ffname2, tnum)
        ns['Removal_'+NAME] = np.array( ns['Removal_'+NAME] )
        # print ns['Removal_'+NAME]

    Removal_Head = Removal_Head_PSD + Removal_Head_noPSD
    Removal_Neck = Removal_Neck_PSD + Removal_Neck_noPSD

    ##
    ## Time development of Tot Ca Concentration
    ##

    TCa_Head = Ca_Head # zeros_like(Ca_Head)
    TCa_Neck = Ca_Neck
    TCa_Dend = Ca_Dend
    CaM       = ['N0C1','N0C2','N1C0','N1C1','N1C2','N2C0','N2C1', 'N2C2', 'CBCa']
    CaMCharge = [     1,     2,     1,     2,     3,     2,     3,      4,      1]
    for (T_CaM, C) in zip(CaM, CaMCharge):
        ffname = p.DIRECTORY + "/"+ T_CaM + "_%03i.npy"
        TCa_Head = TCa_Head + MolTimeDev(ffname, tnum, Head_ID) * C
        TCa_Neck = TCa_Neck + MolTimeDev(ffname, tnum, Neck_ID) * C
        TCa_Dend = TCa_Dend + MolTimeDev(ffname, tnum, Dend_ID) * C
 
    ##
    ## Flux (From_Head) :
    ## dCa = -Flux -Uptake + NMDACa
    ## Flux = -dCa -Uptake + NMDACa
    ##

    tmp       = np.concatenate(( np.array([0]), TCa_Head[0:-1] ))
    dCa_Head  = TCa_Head - tmp 
    tmp       = np.concatenate(( np.array([0]), TCa_Neck[0:-1] ))
    dCa_Neck  = TCa_Neck - tmp 
    tmp       = np.concatenate(( np.array([0]), TCa_Dend[0:-1] ))
    dCa_Dend  = TCa_Dend - tmp 

    Flux_Head = -dCa_Head - Removal_Head + NMDA_Head_Ca
    Flux_Neck1 = - ( dCa_Head + dCa_Neck ) - ( Removal_Head + Removal_Neck ) + ( NMDA_Head_Ca + NMDA_Neck_Ca )
    Flux_Neck2 = dCa_Dend + Removal_Dend

    ##
    ## Plot of time developments (Head/Neck/Dend) 
    ##

    NMDA_Ca = (NMDA_Head_Ca + NMDA_Neck_Ca)
    NMDA_Ca = np.sum(NMDA_Ca)
    Flux_Head = np.sum(Flux_Head)
    Flux_Neck1 = np.sum(Flux_Neck1)


    ##
    ## Save Ca2+, CaM concs
    ##
    ffname1 = p.DIRECTORY + 'Ca_Head.npy'
    ffname2 = p.DIRECTORY + 'Ca_Neck.npy'
    ffname3 = p.DIRECTORY + 'Ca_Dend.npy'
    ffname4 = p.DIRECTORY + 'CaM_Head.npy'
    ffname5 = p.DIRECTORY + 'CaM_Neck.npy'
    ffname6 = p.DIRECTORY + 'CaM_Dend.npy'

    np.save( ffname1, Ca_Head )
    np.save( ffname2, Ca_Neck )
    np.save( ffname3, Ca_Dend )
    np.save( ffname4, CaM_Head )
    np.save( ffname5, CaM_Neck )
    np.save( ffname6, CaM_Dend )


    return NMDA_Ca, Flux_Head,  Flux_Neck1

    ##
    ##
    ##
class Post1_flux_conc():
    def __init__(self, TargDIR):
    ##
    ##
    ##
        p       = parameters( 0 )
        NUM     = p.IDs
    ##
    ## Specify target
    ##
        NMDA_flux = []
        Head_flux = []
        Neck_flux = []

        for i in NUM:
            p   = parameters( i )
            i_str = '%s' % i
            p.DIRECTORY  =  TargDIR + '/' + i_str + '/'
            NMDA_Ca,  Flux_Head,  Flux_Neck = CaFlux( p )
            NMDA_flux.append(NMDA_Ca)
            Head_flux.append(Flux_Head)
            Neck_flux.append(Flux_Neck)

        Neck_NMDAR_ratio = np.array(Neck_flux) / np.array(NMDA_flux)
        Neck_NMDAR_ratio = Neck_NMDAR_ratio.tolist()

        ffname1 = TargDIR + '/NUM.npy'
        ffname2 = TargDIR + '/NMDA_flux.npy'
        ffname3 = TargDIR + '/Head_flux.npy'
        ffname4 = TargDIR + '/Neck_flux.npy'
        ffname5 = TargDIR + "/Neck_NMDAR_ratio.npy"

        np.save( ffname1, NUM )
        np.save( ffname2, NMDA_flux )
        np.save( ffname3, Head_flux )
        np.save( ffname4, Neck_flux )
        np.save( ffname5, Neck_NMDAR_ratio )

        print 'NUM      : ' , NUM
        print 'NMDA_flux: ' , NMDA_flux
        print 'Head_flux: ' , Head_flux
        print 'Neck_flux: ' , Neck_flux

## import sys
## sys.path.append('./POST_PROCESSES')
## from   Post1_flux_conc  import *
## from   Post3_integCa  import *

## Post1_flux_conc('Data/Data4')
## Post3_integCa('Data/Data4')
