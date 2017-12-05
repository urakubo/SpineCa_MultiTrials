# coding: utf-8

import steps.model as smod
import steps.utilities.meshio as smeshio
import steps.geom as sgeom
import math
import numpy as np
import sys
import subprocess
import os
import pickle

####
####
class parameters():
####
####
    def __init__(self, ID):
####
####


        ###
        ### Define geometry data files (ID specific)
        ###
        ID_str = '%s' % ID
        self.MorphID = ID_str +'/170626_2136'

        if ID == 0:
            # self.IDs = [ \
            #    1, 3, 5, 7, 9, 11, 13, 15, 19, 23, 27, 29, 31, 33, \
            #        2, 4, 6, 8, 10, 12, 14, 16, 22, 25, 28, 30, 32, 34]
            #
            # Postprocess
            #
            self.IDs = [ \
                1, 3, 5, 7, 9, 11, 13, 15, 19, 23, 27, 29, 31, 33, \
                    2, 4, 6, 10, 12, 14, 16, 22, 25, 28, 30, 32, 34]
            return

        self.DIRECTORY       = 'Data/Data1/'+ID_str +'/'
        self.GeoFileName     = './Morph/' + self.MorphID+ '_nkf-Lu_STRI3.inp'       
        self.HeadIDFileName  = './Morph/' + self.MorphID+ 'Head_ID.dat'       
        self.NeckIDFileName  = './Morph/' + self.MorphID+ 'Neck_ID.dat'       
        self.DendIDFileName  = './Morph/' + self.MorphID+ 'Dendrite_ID.dat'
        self.INT_DOMAIN_FILE = './Morph/' + self.MorphID+ 'INT_DOMAIN.pickle'

        ###
        ### Examine whether SER domain exists, SER_NAME (ID specific)
        ###
        try:
            with open(self.INT_DOMAIN_FILE, mode='r') as f:
                self.INT_DOMAIN = pickle.load(f)
        except IOError:
            print '"%s" cannot be opened. Check SER/[] by main_MorphSpec_Preprocess.' % self.INT_DOMAIN_FILE
        ###
        ###
        ###

        self.mol = ['Ca']
        self.len_mol = len( self.mol )
        self.scale = 1e-6 * 0.01333
        self.CYT_NAME = 'MeshTetra1'
        self.SER_NAME = 'MeshTetra2'
        self.PSD_NAME = 'MeshTri0'
        # self.DT     = 0.004     # Sampling time-step
        self.DT     = 0.01     # Sampling time-step
        self.INT    = 1.5       # Sim endtime
        self.TIME   = "1500 ms"

        # self.XOff = -0.44638338675 +  5.54488229579
        # self.YOff = -2.19437245075 + 8.27877164341
        # self.ZOff = +2.7085852553 - 2.12963016369
        self.XOff = 3.25907791262
        self.YOff = 2.98384669128
        self.ZOff = 0.0810937530688

        self.DCa          = 223 * 1e-12 ## Sould be specified in SI
        self.DCB          = 20 * 1e-12 ## Sould be specified in SI
        self.DCaM         = 10 * 1e-12 ## Sould be specified in SI
        self.KCST_b       = 2.0e-9  # The reaction constant, production

        # Reactions
        self.kon_CB   = 1e6 * 75
        self.kof_CB   = 1   * 29.5
        self.kon_TN   = 1e6 * 770
        self.kof_TN   = 1   * 160000
        self.kon_RN   = 1e6 * 32000
        self.kof_RN   = 1   * 22000
        self.kon_TC   = 1e6 * 84
        self.kof_TC   = 1   * 2600
        self.kon_RC   = 1e6 * 25
        self.kof_RC   = 1   * 6.5



        # Cytosolic molecules
        self.CYT_MOL_NAME  = ['Ca','N0C0','N0C1','N0C2','N1C0','N1C1','N1C2',
                         'N2C0','N2C1','N2C2','CB','CBCa']

        # Surface molecules
        self.SUR_MOL_NAME  = ['PA','PA_Ca','Leak','NR','NR_Glu','NR_O']

        # Diffusion molecules
        self.D_MOL_NAME    = ['D_Ca','D_N0C0','D_N0C1','D_N0C2','D_N1C0','D_N1C1','D_N1C2',
                              'D_N2C0','D_N2C1','D_N2C2','D_CB','D_CBCa']
        # Surface patches 
        self.SURF_DOMAIN_NAMES = ['Dend', 'Neck_noPSD', 'Head_noPSD','Neck_PSD','Head_PSD']
        self.PSD_DOMAIN_NAMES  = ['Neck_PSD','Head_PSD']

####
####
def gen_model(x, flux):
####
####
    print "Create sim"
    mdl  = smod.Model()

    # Defining chemical compartments
    vsys = smod.Volsys('vsys_Cyt',mdl)
    ssys = smod.Surfsys('ssys',mdl)

    ##
    ## Cytosolic molecule species
    ##
    CYT_MOL_NUM   = len(x.CYT_MOL_NAME)
    CYT_MOL       = []
    for i in range(CYT_MOL_NUM):
        print i, x.CYT_MOL_NAME[i]
        ns = globals()
        ns[x.CYT_MOL_NAME[i]] = smod.Spec(x.CYT_MOL_NAME[i], mdl)
        # CYT_MOL.append( smod.Spec(CYT_MOL_NAME[i], mdl) )

    ##
    ## Surface molecule species
    ##
    SUR_MOL_NUM   = len(x.SUR_MOL_NAME)
    SUR_MOL       = []
    for i in range(SUR_MOL_NUM):
        print i, x.SUR_MOL_NAME[i]
        ns = globals()
        ns[x.SUR_MOL_NAME[i]] = smod.Spec(x.SUR_MOL_NAME[i], mdl)

    ##
    ## Cytosolic molecule reactions
    ##
    smod.Reac('on_CB', vsys, lhs=[Ca,CB], rhs=[CBCa], kcst = x.kon_CB)
    smod.Reac('of_CB', vsys, lhs=[CBCa], rhs=[Ca,CB], kcst = x.kof_CB)

    smod.Reac('on_TN1', vsys, lhs=[Ca,N0C0], rhs=[N1C0], kcst = 2*x.kon_TN)
    smod.Reac('of_TN1', vsys, lhs=[N1C0], rhs=[Ca,N0C0], kcst = x.kof_TN)
    smod.Reac('on_RN1', vsys, lhs=[Ca,N1C0], rhs=[N2C0], kcst = x.kon_RN)
    smod.Reac('of_RN1', vsys, lhs=[N2C0], rhs=[Ca,N1C0], kcst = 2*x.kof_RN)

    smod.Reac('on_TN2', vsys, lhs=[Ca,N0C1], rhs=[N1C1], kcst = 2*x.kon_TN)
    smod.Reac('of_TN2', vsys, lhs=[N1C1], rhs=[Ca,N0C1], kcst = x.kof_TN)
    smod.Reac('on_RN2', vsys, lhs=[Ca,N1C1], rhs=[N2C1], kcst = x.kon_RN)
    smod.Reac('of_RN2', vsys, lhs=[N2C1], rhs=[Ca,N1C1], kcst = 2*x.kof_RN)

    smod.Reac('on_TN3', vsys, lhs=[Ca,N0C2], rhs=[N1C2], kcst = 2*x.kon_TN)
    smod.Reac('of_TN3', vsys, lhs=[N1C2], rhs=[Ca,N0C2], kcst = x.kof_TN)
    smod.Reac('on_RN3', vsys, lhs=[Ca,N1C2], rhs=[N2C2], kcst = x.kon_RN)
    smod.Reac('of_RN3', vsys, lhs=[N2C2], rhs=[Ca,N1C2], kcst = 2*x.kof_RN)

    smod.Reac('on_TC1', vsys, lhs=[Ca,N0C0], rhs=[N0C1], kcst = 2*x.kon_TC)
    smod.Reac('of_TC1', vsys, lhs=[N0C1], rhs=[Ca,N0C0], kcst = x.kof_TC)
    smod.Reac('on_RC1', vsys, lhs=[Ca,N0C1], rhs=[N0C2], kcst = x.kon_RC)
    smod.Reac('of_RC1', vsys, lhs=[N0C2], rhs=[Ca,N0C1], kcst = 2*x.kof_RC)

    smod.Reac('on_TC2', vsys, lhs=[Ca,N1C0], rhs=[N1C1], kcst = 2*x.kon_TC)
    smod.Reac('of_TC2', vsys, lhs=[N1C1], rhs=[Ca,N1C0], kcst = x.kof_TC)
    smod.Reac('on_RC2', vsys, lhs=[Ca,N1C1], rhs=[N1C2], kcst = x.kon_RC)
    smod.Reac('of_RC2', vsys, lhs=[N1C2], rhs=[Ca,N1C1], kcst = 2*x.kof_RC)

    smod.Reac('on_TC3', vsys, lhs=[Ca,N2C0], rhs=[N2C1], kcst = 2*x.kon_TC)
    smod.Reac('of_TC3', vsys, lhs=[N2C1], rhs=[Ca,N2C0], kcst = x.kof_TC)
    smod.Reac('on_RC3', vsys, lhs=[Ca,N2C1], rhs=[N2C2], kcst = x.kon_RC)
    smod.Reac('of_RC3', vsys, lhs=[N2C2], rhs=[Ca,N2C1], kcst = 2*x.kof_RC)

    ##
    ## Ca pump and VGCC
    ##
   
    smod.SReac('PA_To_PA_Ca', ssys, slhs=[PA],   ilhs=[Ca]   , srhs=[PA_Ca], kcst = 150*1e6)
    smod.SReac('PA_Ca_To_PA', ssys, slhs=[PA_Ca]             , srhs=[PA]   , kcst = 12)
    smod.SReac('Leak_To_Leak_Ca', ssys, slhs=[Leak], srhs=[Leak] , irhs=[Ca]   , kcst = 0.015)
    smod.SReac('NR_Glu_To_NR',    ssys, slhs=[NR_Glu], srhs=[NR]     , kcst = 50)
    smod.SReac('NR_Glu_To_NR_O',  ssys, slhs=[NR_Glu], srhs=[NR_O]   , kcst = 200)
    smod.SReac('NR_O_To_NR_Glu',  ssys, slhs=[NR_O]  , srhs=[NR_Glu] , kcst = 50)
    smod.SReac('NR_O_To_NR_O_Ca', ssys, slhs=[NR_O]  , srhs=[NR_O]   , irhs=[Ca], kcst = flux)

    ##
    ## Species on left hand side:  Surface (slhs), Outer comp (olhs), Inner comp (ilhs)
    ## Species on right hand side: Surface (srhs), Outer comp (orhs), Inner comp (irhs)
    ##

    # Diffusion rules
    smod.Diff('D_Ca',   vsys ,Ca, x.DCa)
    smod.Diff('D_N0C0', vsys ,N0C0, x.DCaM)
    smod.Diff('D_N0C1', vsys ,N0C1, x.DCaM)
    smod.Diff('D_N0C2', vsys ,N0C2, x.DCaM)
    smod.Diff('D_N1C0', vsys ,N1C0, x.DCaM)
    smod.Diff('D_N1C1', vsys ,N1C1, x.DCaM)
    smod.Diff('D_N1C2', vsys ,N1C2, x.DCaM)
    smod.Diff('D_N2C0', vsys ,N2C0, x.DCaM)
    smod.Diff('D_N2C1', vsys ,N2C1, x.DCaM)
    smod.Diff('D_N2C2', vsys ,N2C2, x.DCaM)
    smod.Diff('D_CB',   vsys ,CB,   x.DCB)
    smod.Diff('D_CBCa', vsys ,CBCa, x.DCB)

    return mdl

####
####
class gen_geom():
####
####
    def __init__(self, x):

        if x.INT_DOMAIN == []:
            DOMAIN_NAMES = [ x.CYT_NAME, x.PSD_NAME ]
        else:
            DOMAIN_NAMES = [ x.CYT_NAME, x.SER_NAME, x.PSD_NAME ]
            
        print x.GeoFileName
        print x.scale
        print DOMAIN_NAMES
        print 'importAbaqus ... \n'
        mesh, nodeproxy, tetproxy, triproxy = smeshio.importAbaqus(x.GeoFileName, x.scale, DOMAIN_NAMES)
        print 'OK\n'
        print 'Load Head, Neck, and Dend IDs\n'
        Head_ID = np.loadtxt(x.HeadIDFileName) - 1
        Neck_ID = np.loadtxt(x.NeckIDFileName) - 1
        Dend_ID = np.loadtxt(x.DendIDFileName) - 1
        Head_ID = Head_ID.astype(np.int)
        Neck_ID = Neck_ID.astype(np.int)
        Dend_ID = Dend_ID.astype(np.int)
        print 'OK\n'

        ntets      = mesh.ntets
        tet_groups = tetproxy.blocksToGroups()

        Cyt_tets   = tet_groups[x.CYT_NAME]
        nCyt_tets  = len(Cyt_tets)

        if x.INT_DOMAIN != []:
            SER_tets   = tet_groups[x.SER_NAME]
            nSER_tets  = len(SER_tets)
            print 'nGol_tets: ', nSER_tets

        tri_groups = triproxy.blocksToGroups()
        tri_PSD    = tri_groups[x.PSD_NAME]
        tri_nPSD   =  len(tri_PSD)

    ###
    ### Volume system
    ### Annotation of geometry
    ###
        Cytosol = sgeom.TmComp('Cytosol', mesh, Cyt_tets)
        Cytosol.addVolsys('vsys_Cyt')

        #    Golgi = sgeom.TmComp('Golgi', mesh, Gol_tets)
        #    Golgi.addVolsys('vsys_Gol')

    ####
    #### Surface & PSD mesh preparation
    ####

    ## Obtain IDs of neighboring triangles from cytosolic tetrahedrons
        tri_S      = set( mesh.getSurfTris() )
        if x.INT_DOMAIN != []:
            tri_SER   = set()
            for i in SER_tets:
                tri = mesh.getTetTriNeighb(i)
                tri_SER.add(tri[0])
                tri_SER.add(tri[1])
                tri_SER.add(tri[2])
                tri_SER.add(tri[3])
            tri_S.difference(tri_SER)

        ###
        ###
        tri_S_PSD  = tri_S.intersection(tri_PSD)

        tri_Head   = set()
        for i in Head_ID:
            tri = mesh.getTetTriNeighb(i)
            tri_Head.add(tri[0])
            tri_Head.add(tri[1])
            tri_Head.add(tri[2])
            tri_Head.add(tri[3])

        tri_Neck   = set()
        for i in Neck_ID:
            tri = mesh.getTetTriNeighb(i)
            tri_Neck.add(tri[0])
            tri_Neck.add(tri[1])
            tri_Neck.add(tri[2])
            tri_Neck.add(tri[3])

        tri_Dend   = set()
        for i in Dend_ID:
            tri = mesh.getTetTriNeighb(i)
            tri_Dend.add(tri[0])
            tri_Dend.add(tri[1])
            tri_Dend.add(tri[2])
            tri_Dend.add(tri[3])

        tri_S_Head = tri_S.intersection(tri_Head)
        tri_S_Neck = tri_S.intersection(tri_Neck)
        tri_S_Dend = tri_S.intersection(tri_Dend)

        tri_S_Head_PSD   = tri_S_Head.intersection(tri_PSD)
        tri_S_Neck_PSD   = tri_S_Neck.intersection(tri_PSD)
        tri_S_Head_noPSD = tri_S_Head.difference(tri_PSD)
        tri_S_Neck_noPSD = tri_S_Neck.difference(tri_PSD)

    ##
    ## Make list from set
    ##

        tri_S            = list( tri_S )
        tri_S_Dend       = list( tri_S_Dend )
        tri_S_Head_PSD   = list( tri_S_Head_PSD )
        tri_S_Neck_PSD   = list( tri_S_Neck_PSD )
        tri_S_Head_noPSD = list( tri_S_Head_noPSD )
        tri_S_Neck_noPSD = list( tri_S_Neck_noPSD )

    ##
    ## Surface system
    ## Annotation of geometry
    ##
        p_Dend       = sgeom.TmPatch('Dend',      mesh, tri_S_Dend       , icomp = Cytosol)
        p_Head_PSD   = sgeom.TmPatch('Head_PSD',  mesh, tri_S_Head_PSD   , icomp = Cytosol)
        p_Neck_PSD   = sgeom.TmPatch('Neck_PSD',  mesh, tri_S_Neck_PSD   , icomp = Cytosol)
        p_Head_noPSD = sgeom.TmPatch('Head_noPSD',mesh, tri_S_Head_noPSD , icomp = Cytosol)
        p_Neck_noPSD = sgeom.TmPatch('Neck_noPSD',mesh, tri_S_Neck_noPSD , icomp = Cytosol)

        p_Dend.addSurfsys('ssys')
        p_Head_PSD.addSurfsys('ssys')
        p_Neck_PSD.addSurfsys('ssys')
        p_Head_noPSD.addSurfsys('ssys')
        p_Neck_noPSD.addSurfsys('ssys')
 
        self.mesh       = mesh
        self.ntets      = ntets
        self.Cyt_tets   = Cyt_tets
        self.nCyt_tets  = nCyt_tets
        self.Cytosol    = Cytosol

        self.tri_S            = tri_S
        self.tri_S_Dend       = tri_S_Dend
        self.tri_S_Head_PSD   = tri_S_Head_PSD
        self.tri_S_Neck_PSD   = tri_S_Neck_PSD
        self.tri_S_Head_noPSD = tri_S_Head_noPSD
        self.tri_S_Neck_noPSD = tri_S_Neck_noPSD

        self.p_Dend       = p_Dend
        self.p_Head_PSD   = p_Head_PSD
        self.p_Neck_PSD   = p_Neck_PSD
        self.p_Head_noPSD = p_Head_noPSD
        self.p_Neck_noPSD = p_Neck_noPSD

    ##
    ## SER surf mesh
    ##
        if x.INT_DOMAIN != []:
 
            tris_Cyt   = set()
            for i in range(nCyt_tets):
                tris   = mesh.getTetTriNeighb(Cyt_tets[i])
                tris_Cyt.add(tris[0])
                tris_Cyt.add(tris[1])
                tris_Cyt.add(tris[2])
                tris_Cyt.add(tris[3])
            tris_SER   = set()
            for i in range(nSER_tets):
                tris   = mesh.getTetTriNeighb(SER_tets[i])
                tris_SER.add(tris[0])
                tris_SER.add(tris[1])
                tris_SER.add(tris[2])
                tris_SER.add(tris[3])
            tris_S_SER      = tris_Cyt.intersection(tris_SER) 
            self.tris_S_SER =  list(tris_S_SER)
        else:
            self.tris_S_SER = []


    ##
    ## Obtain IDs for fluxes
    ##

        mesh.addROI('Head_ROI', sgeom.ELEM_TET, Head_ID)
        mesh.addROI('Neck_ROI', sgeom.ELEM_TET, Neck_ID)
        mesh.addROI('Dend_ROI', sgeom.ELEM_TET, Dend_ID)

        print 'tri_S  : ', len( tri_S )
        print 'tri_S_Head_PSD : ', len( tri_S_Head_PSD )
        print 'tri_S_Neck_PSD : ', len( tri_S_Neck_PSD )
        print 'tri_S_Head_noPSD : ', len( tri_S_Head_noPSD )
        print 'tri_S_Neck_noPSD : ', len( tri_S_Neck_noPSD )
        print 'tri_S_Dend : ', len( tri_S_Dend )
        tmp = len( tri_S_Head_PSD  + tri_S_Neck_PSD + tri_S_Dend + tri_S_Neck_noPSD + tri_S_Head_noPSD )
        print 'Head_PSD + Neck_PSD + Dend + Neck_noPSD + Head_noPSD: ', tmp

#        print  mesh.getSurfSys()
        print 'Gend_geom has been finished!\n'

#        print  p_Dend.getAllPatches()



