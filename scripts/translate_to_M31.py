# This script reads the unique set of complete triplets and translates the data
# to the M31-centric frame.

import numpy as np
import sys
sys.path.append('/afs/slac.stanford.edu/u/ki/yymao/scripts')
from readHlist import readHlist

path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/data_files/all_unique_triplets.txt'
triplets_data = readHlist(path)
header = '# idMW idM31 idM33 MW_Mvir M31_Mvir M33_Mvir MW_x MW_y MW_z M31_x M31_y M31_z M33_x M33_y M_33_z\
 MW_vx MW_vy MW_vz M31_vx M31_vy M31_vz M33_vx M33_vy M33_vz\n'

processed_triplets = np.ndarray(triplets_data.shape[0], 6)

 
