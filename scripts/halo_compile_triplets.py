# This script gathers all the complete triplets dispersed in the 50 Consuelo snapshots 
# and writes one file containing all of them. A complete triplet is one with an M33.

import numpy as np
import sys
sys.path.append('/afs/slac.stanford.edu/u/ki/yymao/scripts')
from readHlist import readHlist

header = '# idMW idM31 idM33 Mvir x y z vx vy vz\n'

path = '/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Consuelo/LG_project/'
prefix = 'MW_halos_'
suffix = '_z0p000_vir_bgc2.txt'

Triplet_file = open('MW_all_Triplets_z0p000_data.txt','a+')
Triplet_file.write(header)
for i in range(4001,4051):
    halo_file_name = path+prefix+str(i)+suffix
    halo_data = readHlist(halo_file_name)
    print halo_file_name
    for data_line in halo_data:
       # print (data_line)
        if data_line['idM33']!=-1.0:
            Triplet_file.write("%d %d %d %f %f %f %f %f %f %f\n"%(data_line[0], data_line[1], data_line[2], data_line[3], data_line[4], data_line[5], data_line[6], data_line[7], data_line[8], data_line[9]))

Triplet_file.close()
