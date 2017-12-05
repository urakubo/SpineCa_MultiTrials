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

def plot3comp2(ax, tpnts, Ca_Head, Ca_Neck, Ca_Dend):
    ms = 3
    ax.plot(tpnts, Ca_Head,'o', markersize=ms, markerfacecolor=b     , markeredgecolor=b     , label='Head')
    ax.plot(tpnts, Ca_Neck,'o', markersize=ms, markerfacecolor=lb    , markeredgecolor=lb    , label='Neck')
    ax.plot(tpnts, Ca_Dend,'o', markersize=ms, markerfacecolor='0.5' , markeredgecolor='0.5' , label='Dend')

def xyrange(ax,xr,yr,xtitle,ytitle):
    ax.set_xlim(xr)
    ax.set_ylim(yr)    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.get_xaxis().set_tick_params(which='both',direction='out')
    ax.get_yaxis().set_tick_params(which='both',direction='out')

def xxrange(ax,xr,xtitle,ytitle):
    ax.set_xlim(xr)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    ax.get_xaxis().set_tick_params(which='both',direction='out')
    ax.get_yaxis().set_tick_params(which='both',direction='out')

def mov_av(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


##
##
##
if __name__ == '__main__':
##
##
##
    TargDIR_Tot = [ 'Data1' , 'Data2', 'Data3', 'Data4' , 'Data5', 'Data6' ]
    argvs = sys.argv
    ##
    ##
    NUM = int(argvs[1])
    p      = parameters( 0 )
    ID     = p.IDs[NUM]
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
    ## Plot of time developments (Head/Neck/Dend) 
    ##
    xtitle = 'Time (s)'
    xr     = [tpnts.min(), tpnts.max()]
    lb     = (0.6, 0.6, 1)
    b      = (0, 0, 1)
    FS     = 11
    fig = plt.figure()
    fig.patch.set_facecolor('w')
    ##
    ax1  = fig.add_subplot(3,3,1)
    ax1.set_title("Free Ca")
    xxrange(ax1, xr, xtitle, 'uM')

    ax2  = fig.add_subplot(3,3,2)
    ax2.set_title("Free Ca")
    xxrange(ax2, xr, xtitle, 'uM')
    ##
    ax3  = fig.add_subplot(3,3,7)
    ax3.set_title("Free CaM")
    xxrange(ax3, xr, xtitle, 'uM')

    ax4  = fig.add_subplot(3,3,8)
    ax4.set_title("Free CaM")
    xxrange(ax4, xr, xtitle, 'uM')

    ##
    ##
    ##
    PeakCaHead  = []
    PeakCaNeck  = []
    PeakCaDend  = []
    PeakCaMHead = []
    PeakCaMNeck = []
    PeakCaMDend = []
    PeakSmoothCaHead  = []
    PeakSmoothCaNeck  = []
    PeakSmoothCaDend  = []
    PeakSmoothCaMHead = []
    PeakSmoothCaMNeck = []
    PeakSmoothCaMDend = []
    ##
    ##
    for TargDIR in TargDIR_Tot:
    ##
    ##

        p.DIRECTORY = 'Data/' + TargDIR + '/' + ID_str + '/'

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

        Ca_Head_C = Ca_Head/Head_Vol[NUM]/(6.02e20)
        Ca_Neck_C = Ca_Neck/Neck_Vol[NUM]/(6.02e20)
        Ca_Dend_C = Ca_Dend/Dend_Vol[NUM]/(6.02e20)

        CaM_Head_C = CaM_Head/Head_Vol[NUM]/(6.02e20)
        CaM_Neck_C = CaM_Neck/Neck_Vol[NUM]/(6.02e20)
        CaM_Dend_C = CaM_Dend/Dend_Vol[NUM]/(6.02e20)

        PeakCaHead.append( np.max( Ca_Head_C ) )
        PeakCaNeck.append( np.max( Ca_Neck_C ) )
        PeakCaDend.append( np.max( Ca_Dend_C ) )
        PeakCaMHead.append( np.max( CaM_Head_C ) )
        PeakCaMNeck.append( np.max( CaM_Neck_C ) )
        PeakCaMDend.append( np.max( CaM_Dend_C ) )

        plot3comp(ax1, tpnts, Ca_Head_C, Ca_Neck_C, Ca_Dend_C)
        plot3comp(ax3, tpnts, CaM_Head_C, CaM_Neck_C, CaM_Dend_C)

        num = 5
        Ca_Head_C  = mov_av(Ca_Head_C, num )
        Ca_Neck_C  = mov_av(Ca_Neck_C, num )
        Ca_Dend_C  = mov_av(Ca_Dend_C, num )
        CaM_Head_C = mov_av(CaM_Head_C, num )
        CaM_Neck_C = mov_av(CaM_Neck_C, num )
        CaM_Dend_C = mov_av(CaM_Dend_C, num )
        ##
        PeakSmoothCaHead.append(  np.max( Ca_Head_C ) )
        PeakSmoothCaNeck.append(  np.max( Ca_Neck_C ) )
        PeakSmoothCaDend.append(  np.max( Ca_Dend_C ) )
        PeakSmoothCaMHead.append( np.max( CaM_Head_C ) )
        PeakSmoothCaMNeck.append( np.max( CaM_Neck_C ) )
        PeakSmoothCaMDend.append( np.max( CaM_Dend_C ) )
        ##

        plot3comp(ax2, tpnts, Ca_Head_C, Ca_Neck_C, Ca_Dend_C)
        plot3comp(ax4, tpnts, CaM_Head_C, CaM_Neck_C, CaM_Dend_C)

        ##

    ones = list( np.ones(len(TargDIR_Tot)) ) 
    twos = list( np.ones(len(TargDIR_Tot))*2 ) 
        ##
    ax1  = fig.add_subplot(3,3,3)
    ax1.set_title("Free Ca")
    xxrange(ax1, [0.5, 2.5], xtitle, 'uM')
    plot3comp2(ax1, ones , PeakCaHead, PeakCaNeck, PeakCaDend)
    plot3comp2(ax1, twos , PeakSmoothCaHead, PeakSmoothCaNeck, PeakSmoothCaDend)
    ax1.set_xticks([1,2])
    ax1.set_xticklabels(['Org','Smooth'])
    

    ax2  = fig.add_subplot(3,3,9)
    ax2.set_title("Free CaM")
    xxrange(ax2, [0.5, 2.5], xtitle, 'uM')
    plot3comp2(ax2, ones , PeakCaMHead, PeakCaMNeck, PeakCaMDend)
    plot3comp2(ax2, twos , PeakSmoothCaMHead, PeakSmoothCaMNeck, PeakSmoothCaMDend)
    ax2.set_xticks([1,2])
    ax2.set_xticklabels(['Org','Smooth'])

    plt.show()

##################
##################
##################
##################
##################



