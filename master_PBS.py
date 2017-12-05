# coding: utf-8

###
### Usage:
### mpirun -n 2 python main170616_mpi.py


from sub171111 import *

import os


p    = parameters( 0 )


for j in range(1,4):
    print j
    for i in p.IDs:
        tmp = "echo 'source ~/.bashrc; cd " + (os.getcwd())
        cmd = tmp + ";python worker_PBS.py {0} {1}' | qsub -lnodes=maui46 -o PBS_OUTPUT -e PBS_OUTPUT".format(j,i)
        print cmd
        os.system(cmd)

for j in range(4,7):
    print j
    for i in p.IDs:
        tmp = "echo 'source ~/.bashrc; cd " + (os.getcwd())
        cmd = tmp + ";python worker_PBS.py {0} {1}' | qsub -lnodes=maui47 -o PBS_OUTPUT -e PBS_OUTPUT".format(j,i)
        print cmd
        os.system(cmd)



#
# qdelall
# Outout dir of STDIN should be specified by :
# -o ~/LaxmiSpine/170618_Str/PBS_OUTPUT
#
