# This script gathers all the complete triplets dispersed in the 50 Consuelo snapshots 
# and writes one file containing all of them. A complete triplet is one with an M33.

import numpy as np
import sys
sys.path.append('/afs/slac.stanford.edu/u/ki/yymao/scripts')
from readHlist import readHlist

# idMW idM31 idM33 Mvir x y z vx vy vz
header = '# idMW idM31 idM33 MW_Mvir M31_Mvir M33_Mvir MW_x MW_y MW_z M31_x M31_y M31_z M33_x M33_y M_33_z\
 MW_vx MW_vy MW_vz M31_vx M31_vy M31_vz M33_vx M33_vy M33_vz\n'
path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/data_files/'
M33_file = 'M33_all_z0p000_data.txt'
MW_trip_file = 'MW_all_Triplets_z0p000_data.txt'
M31_trip_file = 'M31_all_Triplets_z0p000_data.txt'

M33_data = readHlist(path+M33_file)
MW_data = readHlist(path+MW_trip_file)
M31_data = readHlist(path+M31_trip_file)


triplet_file = open('complete_triplets.txt', 'w')
triplet_file.write(header)
for MW_line in MW_data:
    MW_id = MW_line['idMW']
    M31_id = MW_line['idM31']
    M33_id = MW_line['idM33']
    M31_line = M31_data[M31_data['idM31']==M31_id][0]
    M33_line = M33_data[M33_data['idM33']==M33_id][0]
    triplet_file.write('%d %d %d %f %f %f '%(MW_line['idMW'], M31_line['idM31'], M33_line['idM33'], \
MW_line['Mvir'], M31_line['Mvir'], M33_line['Mvir']))
    triplet_file.write('%f %f %f %f %f %f %f %f %f '%(MW_line['x'], MW_line['y'], MW_line['z'], \
M31_line['x'], M31_line['y'], M31_line['z'], M33_line['x'], M33_line['y'], M33_line['z']))
    triplet_file.write('%f %f %f %f %f %f %f %f %f\n'%(MW_line['vx'], MW_line['vy'], MW_line['vz'], \
M31_line['vx'], M31_line['vy'], M31_line['vz'], M33_line['vx'], M33_line['vy'], M33_line['vz']))
