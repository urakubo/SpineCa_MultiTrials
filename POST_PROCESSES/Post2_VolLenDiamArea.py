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


##
## python -i main_plot.py
##
   
def getVol(targ_mesh, ID):
    Vol = 0.0
    for i in ID:
        Vol = Vol + targ_mesh.getTetVol(i)
    return Vol

def VolHist(targ_mesh, xx, Head_ID, Head_Loc, Neck_ID, Neck_Loc):
    Locs   = []
    Vols   = []
    for (i, id) in enumerate(Head_ID):
        Vols.append( targ_mesh.getTetVol(id) )
        Locs.append( Head_Loc[i] )
    for (i, id) in enumerate(Neck_ID):
        Vols.append( targ_mesh.getTetVol(id) )
        Locs.append( Neck_Loc[i] )
    [tmpH, xx] = np.histogram(Locs, xx, weights=Vols)
    return tmpH


def MolHist(Mol, xx, Tj, Head_ID, Head_Loc, Neck_ID, Neck_Loc):
    Hj = np.zeros(( xx.size-1, Tj.size ))
    for (k,j) in enumerate(Tj):
        ffname = p.DIRECTORY[0] + "/" + Mol + "_%03i.npy" % j
        Ca = np.load(ffname)
        locs   = []
        for (i, id) in enumerate(Head_ID):
            tmpCa  = Ca[id]
            tmploc = Head_Loc[i]
            for jj in range(np.int(tmpCa)):
                locs.append(tmploc)
        for (i, id) in enumerate(Neck_ID):
            tmpCa  = Ca[id]
            tmploc = Neck_Loc[i]
            for jj in range(np.int(tmpCa)):
                locs.append(tmploc)
        [tmpH, xx] = np.histogram(locs, xx)
        # print tmpH
        Hj[:,k] = tmpH
    return Hj

def MeanSectionArea( p, g, xx ):
    Targ = './Morph/' + p.MorphID
    Head_ID   = np.loadtxt(Targ+'Head_ID.dat') - 1
    Neck_ID   = np.loadtxt(Targ+'Neck_ID.dat') - 1
    Head_Loc  = np.loadtxt(Targ+'Head_Loc.dat')
    Neck_Loc  = np.loadtxt(Targ+'Neck_Loc.dat')
    Head_ID = Head_ID.astype(np.int)
    Neck_ID = Neck_ID.astype(np.int)

    volumeHist = VolHist(g.mesh, xx, Head_ID, Head_Loc, Neck_ID, Neck_Loc)
    a = np.array( volumeHist / 0.1 ) * np.power(10, 18)
    MeanArea   = a.sum()/np.nonzero(a)[0].size
    return MeanArea

def SpineNeckArea( p, g, xx ):
    Targ = './Morph/' + p.MorphID
    Head_ID   = np.loadtxt(Targ+'Head_ID.dat') - 1
    Neck_ID   = np.loadtxt(Targ+'Neck_ID.dat') - 1
    Head_Loc  = np.loadtxt(Targ+'Head_Loc.dat')
    Neck_Loc  = np.loadtxt(Targ+'Neck_Loc.dat')
    Head_ID = Head_ID.astype(np.int)
    Neck_ID = Neck_ID.astype(np.int)

    volumeHist       = VolHist(g.mesh, xx, Head_ID, Head_Loc, Neck_ID, Neck_Loc)
    a = np.array( volumeHist / 0.1 ) * np.power(10, 18)
    a = a[2:]
    a = np.delete(a,np.where(a==0)[0])
    MinNeckArea = np.min( a )
    return MinNeckArea

##
##
##
class Post2_VolLenDiamArea():
    def __init__(self, TargDIR):
    ##
    ##
    ##
        p   = parameters( 0 )
        NUM = p.IDs
    ##
    ##
        HeadVolume   = []
        NeckVolume   = []
        DendVolume   = []
        SpineLen     = []
        MeanNeckArea = []
        MinNeckArea  = []
        HeadArea     = []
        NeckArea     = []
        PSDArea      = []
        xx = np.arange(-0.1, 2.5, 0.1) # start, stop, step
    ##
    ##
    ##
        for i in NUM:
            p   = parameters( i )
            g   = gen_geom(p)
            Targ = './Morph/' + p.MorphID
            print Targ
        ##
        ## Head volume
        ##
            Head_ID   = np.loadtxt(Targ+'Head_ID.dat') - 1
            Head_ID   = Head_ID.astype(np.int)
            Neck_ID   = np.loadtxt(Targ+'Neck_ID.dat') - 1
            Neck_ID   = Neck_ID.astype(np.int)
            Dend_ID   = np.loadtxt(Targ+'Dendrite_ID.dat') - 1
            Dend_ID   = Dend_ID.astype(np.int)
            tmp = getVol(g.mesh, Head_ID)
            HeadVolume.append( tmp )
            tmp = getVol(g.mesh, Neck_ID)
            NeckVolume.append( tmp )
            tmp = getVol(g.mesh, Dend_ID)
            DendVolume.append( tmp )
        ##
        ## Spinelen, MinNeck, MeanNeck
        ##
            Dist_Area = np.loadtxt(Targ+'Dist_Area.dat')
            SpineLen.append(Dist_Area[-1,0])
            tmp = SpineNeckArea( p, g, xx )
            MinNeckArea.append(tmp)
            tmp = MeanSectionArea( p, g, xx )
            MeanNeckArea.append(tmp)
        ##
        ## SurfArea
        ##
            tmp1 = g.p_Neck_noPSD.getArea()*1e12
            tmp2 = g.p_Neck_PSD.getArea()*1e12
            NeckArea.append( tmp1+tmp2 )
            tmp3 = g.p_Head_noPSD.getArea()*1e12
            tmp4 = g.p_Head_PSD.getArea()*1e12
            HeadArea.append( tmp3+tmp4 )
            PSDArea.append( tmp2+tmp4 )

        ffname1 = TargDIR + "/HeadVol.npy"
        ffname2 = TargDIR + "/NeckVol.npy"
        ffname3 = TargDIR + "/DendVol.npy"
        np.save( ffname1, HeadVolume )
        np.save( ffname2, NeckVolume )
        np.save( ffname3, DendVolume )

        ffname2 = TargDIR + "/SpineLen.npy"
        ffname3 = TargDIR + "/MinNeckArea.npy"
        ffname4 = TargDIR + "/MeanNeckArea.npy"
        np.save( ffname2, SpineLen )
        np.save( ffname3, MinNeckArea )
        np.save( ffname4, MeanNeckArea )

        ffname1 = TargDIR + "/HeadArea.npy"
        ffname2 = TargDIR + "/NeckArea.npy"
        ffname3 = TargDIR + "/PSDArea.npy"
        np.save( ffname1, HeadArea )
        np.save( ffname2, NeckArea )
        np.save( ffname3, PSDArea )


