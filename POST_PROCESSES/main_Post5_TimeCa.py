##
## python -i main170622_plot.py
##
from sub171111 import *
import math
import sys
import numpy as np
from itertools import count


def MolTimeDev(fname, tnum, id):
    id = np.int_(id)
    Mol = []
    for j in range(tnum):
        ffname = fname % j
        tmp = sum( np.load(ffname)[id] )
        Mol.append( tmp.tolist() )
    Mol = np.array( Mol )
    return Mol


if __name__ == '__main__':

    ##
    ## Specify target
    ## 
    argvs = sys.argv
    if   argvs[1] == '1':
        TargDIR = 'Data1'
    elif argvs[1] == '2':
        TargDIR = 'Data2'
    elif argvs[1] == '3':
        TargDIR = 'Data3'
    else:
        print 'No simulation data\n'
        sys.exit()
    ##
    ##
    p   = parameters( 0 )
    NUM = p.IDs
    ##
    ## Specify target
    ## 
    NMDA_flux = []
    Head_flux = []
    Neck_flux = []

    for i in NUM:
        p   = parameters( i )
        i_str = '%s' % i
        print 'NUM: '+i_str
        p.DIRECTORY  =  TargDIR + '/' + i_str + '/'
        Targ        = './Morph/' + p.MorphID
    ##
    ##
        tpnts = np.arange(0.0, p.INT, p.DT)
        tpnts = tpnts - tpnts[50]
        tnum  = len( tpnts )
    ##
    ##
        Head_ID   = np.loadtxt(Targ+'Head_ID.dat') - 1
        Neck_ID   = np.loadtxt(Targ+'Neck_ID.dat') - 1
        Dend_ID   = np.loadtxt(Targ+'Dendrite_ID.dat') - 1
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
    ##
    ##
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

##################
##################
