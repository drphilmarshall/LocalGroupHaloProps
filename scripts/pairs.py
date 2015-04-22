# This script gathers all the complete triplets dispersed in the 50 Consuelo snapshots 
# and writes one file containing all of them. A complete triplet is one with an M33.

import numpy as np
import sys

sys.path.append('/u/ki/yymao/pyscripts')
from helpers.SimulationAnalysis import readHlist

header = '# idMW idM31 MW_Mvir M31_Mvir MW_x MW_y MW_z M31_x M31_y M31_z MW_vx MW_vy MW_vz M31_vx M31_vy M31_vz\n'

path = '/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Consuelo/LG_project/'
MW_prefix = 'MW_halos_'
M31_prefix = 'M31_halos_'
suffix = '_z0p000_vir_bgc2.txt'

pair_file = open('MW_M31_pairs.txt','a+')
pair_file.write(header)
for i in range(4001,4051):
    MW_halo_file_name = path+MW_prefix+str(i)+suffix
    MW_data = readHlist(MW_halo_file_name)
    MW_data = MW_data[MW_data['idM33']==-1]
    M31_halo_file_name = path+M31_prefix+str(i)+suffix
    M31_data = readHlist(M31_halo_file_name)
    M31_data = M31_data[M31_data['idM33']==-1]
    print i
    np.random.shuffle(MW_data)
    for MW_line in MW_data[0:30000]:
       # print (data_line)
        MW_id = MW_line['idMW']
        M31_line = M31_data[M31_data['idMW']==MW_id]
        pair_file.write("%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n"%(MW_line['idMW'], MW_line['idM31'], MW_line['Mvir'], M31_line['Mvir'], MW_line['x'], MW_line['y'], MW_line['z'], M31_line['x'], M31_line['y'], M31_line['z'], MW_line['vx'], MW_line['vy'], MW_line['vz'], M31_line['vx'], M31_line['vy'], M31_line['vz']))


pair_file.close()
