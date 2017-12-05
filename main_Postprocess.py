# coding: utf-8

from sub171111 import *
import math
import numpy as np
import sys
import subprocess
import os
import pickle

sys.path.append('./POST_PROCESSES')
from   Post1_flux_conc import *
from   Post2_VolLenDiamArea import *
from   Post3_integCa import *


##
## make Data storage directory
##
TargDIR_Tot = [ 'Data1' , 'Data2', 'Data3', \
                    'Data4' , 'Data5', 'Data6']

for TargDIR in TargDIR_Tot:
    DIR = 'Data/' +  TargDIR
    Post1_flux_conc(DIR)
    Post2_VolLenDiamArea(DIR)
    Post3_integCa(DIR)
