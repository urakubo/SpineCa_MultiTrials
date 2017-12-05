# coding: utf-8

from sub171111 import *
import math
import numpy as np
import sys
import subprocess
import os
import pickle

####
####

##
## make Data storage directory
##
TargDIR_Tot = [ 'Data1' , 'Data2', 'Data3', \
                    'Data4' , 'Data5', 'Data6']

Exec = 'mkdir '+ (os.getcwd()) + '/Data'
os.system(Exec)
for TargDIR in TargDIR_Tot:
    Exec = 'mkdir '+ (os.getcwd()) + '/Data/' +  TargDIR
    print Exec
    os.system(Exec)
    Exec = 'mkdir '+ (os.getcwd()) + '/Data/' +  TargDIR + '/{1..37}'
    print Exec
    os.system(Exec)

####
####
p = parameters( 0 )
INT_DOMAIN = [];
####
####

##
##
for i in p.IDs:
##
##
    p = parameters( i )

##
##  Change char code
##
    IInput = (os.getcwd()) + '/Morph/' + p.MorphID+ '.inp'    
    IOutput = (os.getcwd()) + '/Morph/' + p.MorphID+ '_nkf-Lu_STRI3.inp'
    Exec = '/usr/bin/nkf -Lu '+ IInput + ' > ' + IOutput
    print 'Input Geo File : ', IInput
    print 'Output Geo File: ', IOutput
    print Exec
    os.system(Exec)

##
##  Examine whether internal domain exists, SER_NAME
##
    ld  = open(IOutput)
    lines = ld.readlines()
    ld.close()
    tmp = 0
    for line in lines:
        if line.find( p.SER_NAME ) >= 0:
            tmp = tmp + 1
    if tmp > 0:
        INT_DOMAIN = p.SER_NAME
    else:
        INT_DOMAIN = []
    print i, INT_DOMAIN

##
##  Save ON/OFF of internal domain
##
    with open(p.INT_DOMAIN_FILE, mode='w') as f:
        pickle.dump(INT_DOMAIN, f)


# with open('Morph/INT_DOMAIN.pickle', mode='r') as f:
# tmp = pickle.load(f)
