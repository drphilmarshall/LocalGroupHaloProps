import numpy as np
import sys
import os
sys.path.append('/u/ki/yymao/pyscripts')
from helpers.SimulationAnalysis import readHlist

# idMW idM31 idM33 Mvir x y z vx vy vz
header = '# idMW idM31 MW_Mvir M31_Mvir MW_x MW_y MW_z M31_x M31_y M31_z MW_vx MW_vy MW_vz M31_vx M31_vy M31_vz\n'


path = '/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Consuelo/LG_project/'
MW_prefix = 'MW_halos_'
M31_prefix = 'M31_halos_'
suffix = '_z0p000_vir_bgc2.txt'

pair_file = open('pairs.txt','a+')
pair_file.write(header)
for i in range(4001,4051):
    MW_file = 
    print M33_file_name
    for data_line in M33_data:
