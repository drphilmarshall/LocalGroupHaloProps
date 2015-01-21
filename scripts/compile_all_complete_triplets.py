# This script gathers all the complete triplets dispersed in the 50 Consuelo snapshots 
# and writes one file containing all of them. A complete triplet is one with an M33.

import numpy as np
import sys
sys.path.append('/afs/slac.stanford.edu/u/ki/yymao/scripts')
from readHlist import readHlist

header = '# idMW idM31 idM33 Mvir x y z vx vy vz\n'

path = '/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Consuelo/LG_project/'
M31_prefix = 'M31_halos_'
M31_suffix = '_z0p000_vir_bgc2.txt'

Triplet_file = open('All_Triplets_z0p000_data.txt','a+')
Triplet_file.write(header)
for i in range(4001,4051):
    M31_file_name = path+M31_prefix+str(i)+M31_suffix
    M31_data = readHlist(M31_file_name)
    print M31_file_name
    for data_line in M31_data:
       # print (data_line)
        if data_line['idM33']!=-1.0:
            Triplet_file.write("%d %d %d %f %f %f %f %f %f %f\n"%(data_line[0], data_line[1], data_line[2], data_line[3], data_line[4], data_line[5], data_line[6], data_line[7], data_line[8], data_line[9]))

Triplet_file.close()
