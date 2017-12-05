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


def plot3comp(ax, tpnts, Ca_Head, Ca_Neck, Ca_Dend):
    ax.plot(tpnts, Ca_Head,'-', color=b     , label='Head')
    ax.plot(tpnts, Ca_Neck,'-', color=lb    , label='Neck')
    ax.plot(tpnts, Ca_Dend,'-', color='0.5' , label='Dend')

    
def xyrange(ax,xr,yr,xtitle,ytitle):
    ax.set_xlim(xr)
    ax.set_ylim(yr)    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)


    ##
    ##
if __name__ == '__main__':
    ##
    ##
    TargDIR_Tot = [ 'Data1' , 'Data2', 'Data3' ]
    argvs = sys.argv
    ##
    ##
    p      = parameters( 0 )
    NUM    = p.IDs
    ID     = int(argvs[1])
    p      = parameters( ID )
    ID_str = '%s' % ID
    Targ   = './Morph/' + p.MorphID
    ##
    ##
    tpnts  = np.arange(0.0, p.INT, p.DT)
    tpnts  = tpnts - tpnts[50]
    tnum   = len( tpnts )
    ##
    ##
    ffname1 = 'Data/' + TargDIR_Tot[0] + "/HeadVol.npy"
    ffname2 = 'Data/' + TargDIR_Tot[0] + "/NeckVol.npy"
    ffname3 = 'Data/' + TargDIR_Tot[0] + "/DendVol.npy"
    Head_Vol = np.load( ffname1 )
    Neck_Vol = np.load( ffname2 )
    Dend_Vol = np.load( ffname3 )
    ##
    ##
    Ca_Head_Tot  = np.zeros( tnum )
    Ca_Neck_Tot  = np.zeros( tnum )
    Ca_Dend_Tot  = np.zeros( tnum )
    CaM_Head_Tot = np.zeros( tnum )
    CaM_Neck_Tot = np.zeros( tnum )
    CaM_Dend_Tot = np.zeros( tnum )
    ##
    ##
    for TargDIR in TargDIR_Tot:
    ##
    ##
        p.DIRECTORY = 'Data/' + TargDIR + '/' + argvs[1] + '/'
    ##
    ## Obtain volumes (in m3)
    ##
        ffname1 = p.DIRECTORY + 'Ca_Head.npy'
        ffname2 = p.DIRECTORY + 'Ca_Neck.npy'
        ffname3 = p.DIRECTORY + 'Ca_Dend.npy'
        ffname4 = p.DIRECTORY + 'CaM_Head.npy'
        ffname5 = p.DIRECTORY + 'CaM_Neck.npy'
        ffname6 = p.DIRECTORY + 'CaM_Dend.npy'

        Ca_Head  = np.load( ffname1 )
        Ca_Neck  = np.load( ffname2 )
        Ca_Dend  = np.load( ffname3 )
        CaM_Head = np.load( ffname4 )
        CaM_Neck = np.load( ffname5 )
        CaM_Dend = np.load( ffname6 )

        Ca_Head_C = Ca_Head/Head_Vol[NUM[ID]]/(6.02e20)
        Ca_Neck_C = Ca_Neck/Neck_Vol[NUM[ID]]/(6.02e20)
        Ca_Dend_C = Ca_Dend/Dend_Vol[NUM[ID]]/(6.02e20)

        CaM_Head_C = CaM_Head/Head_Vol[NUM[ID]]/(6.02e20)
        CaM_Neck_C = CaM_Neck/Neck_Vol[NUM[ID]]/(6.02e20)
        CaM_Dend_C = CaM_Dend/Dend_Vol[NUM[ID]]/(6.02e20)

        Ca_Head_Tot  = Ca_Head_C  + Ca_Head_Tot
        Ca_Neck_Tot  = Ca_Neck_C  + Ca_Neck_Tot
        Ca_Dend_Tot  = Ca_Dend_C  + Ca_Dend_Tot
        CaM_Head_Tot = CaM_Head_C + CaM_Head_Tot
        CaM_Neck_Tot = CaM_Neck_C + CaM_Neck_Tot
        CaM_Dend_Tot = CaM_Dend_C + CaM_Dend_Tot
    ##
    ##
    Ca_Head_Tot  = Ca_Head_Tot / len( TargDIR_Tot )
    Ca_Neck_Tot  = Ca_Neck_Tot / len( TargDIR_Tot )
    Ca_Dend_Tot  = Ca_Dend_Tot / len( TargDIR_Tot )
    CaM_Head_Tot = CaM_Head_Tot / len( TargDIR_Tot )
    CaM_Neck_Tot = CaM_Neck_Tot / len( TargDIR_Tot )
    CaM_Dend_Tot = CaM_Dend_Tot / len( TargDIR_Tot )


    ##
    ## Plot of time developments (Head/Neck/Dend) 
    ##

    fig = plt.figure()
    fig.patch.set_facecolor('w')
    xtitle = 'Time (s)'
    xr     = [tpnts.min(), tpnts.max()]
    lb     = (0.6, 0.6, 1)
    b      = (0, 0, 1)
    FS     = 11

    ymx = np.max( np.concatenate((Ca_Head_Tot, Ca_Neck_Tot, Ca_Dend_Tot)) )
    ax  = fig.add_subplot(3,3,1)
    ax.set_title("Ave Free Ca")
    xyrange(ax, xr,[-0.1*ymx,1.1*ymx], xtitle, 'uM')
    plot3comp(ax, tpnts, Ca_Head_Tot, Ca_Neck_Tot, Ca_Dend_Tot)
    ax.legend(frameon=False, fontsize=FS)

    ax  = fig.add_subplot(3,3,7)
    ax.set_title("Free Ca")
    ymx = np.max( np.concatenate((Ca_Head_C, Ca_Neck_C, Ca_Dend_C)) )
    xyrange(ax, xr,[-0.1*ymx,1.1*ymx], xtitle, 'uM')
    plot3comp(ax, tpnts, Ca_Head_C, Ca_Neck_C, Ca_Dend_C)

 
    ymx = np.max( np.concatenate((CaM_Head_Tot, CaM_Neck_Tot, CaM_Dend_Tot)) )   
    ax  = fig.add_subplot(3,3,3)
    ax.set_title("Ave Ca-bound CaM")
    xyrange(ax, xr,[-0.1*ymx,1.1*ymx], xtitle, 'uM')
    plot3comp(ax, tpnts, CaM_Head_Tot, CaM_Neck_Tot, CaM_Dend_Tot)


    ymx = np.max( np.concatenate((CaM_Head_C, CaM_Neck_C, CaM_Dend_C)) )   
    ax  = fig.add_subplot(3,3,9)
    ax.set_title("Ca-bound CaM")
    xyrange(ax, xr,[-0.1*ymx,1.1*ymx], xtitle, 'uM')
    plot3comp(ax, tpnts, CaM_Head_C, CaM_Neck_C, CaM_Dend_C)

    plt.show()

##################
##################
##################
##################
##################

