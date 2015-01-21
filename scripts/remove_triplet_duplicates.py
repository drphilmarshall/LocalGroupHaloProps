# This script should be applied after all the complete triplets have been 
# gathered into one data file using compile_all_complete_triplets.py.
# It turns out that there are rare cases where a single idMW, idM31, or idM33
# is part of multiple triplets. Physically, this case does not resemble our 
# Local Group at all, so it should be excluded.

def filter_unique(data, field):

    if len(data[field]) != len(np.unique(data[field])):
        field_freqs = itemfreq(data[field])
        freq_dict = dict()
        for pair in field_freqs:
            freq_dict[pair[0]] = pair[1]
        keep_data = np.array([freq_dict[entry]==1 for entry in data[field]])
        filtered_data = data[keep_data]
    else:
        filtered_data = data

    return filtered_data

from scipy.stats import itemfreq
import numpy as np
import sys
sys.path.append('/afs/slac.stanford.edu/u/ki/yymao/scripts')
from readHlist import readHlist

header = '# idMW idM31 idM33 Mvir x y z vx vy vz\n'

all_triplet_path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/data_files/All_Triplets_z0p000_data.txt'

all_triplet_data = readHlist(all_triplet_path)

all_triplet_data = filter_unique(all_triplet_data, 'idMW')
all_triplet_data = filter_unique(all_triplet_data, 'idM31')
all_triplet_data = filter_unique(all_triplet_data, 'idM33')

filtered_file = open('all_unique_triplets.txt', 'w')
filtered_file.write(header)
for data_line in all_triplet_data:
    filtered_file.write("%d %d %d %f %f %f %f %f %f %f\n"%(data_line[0], data_line[1], data_line[2], data_line[3], data_line[4], data_line[5], data_line[6], data_line[7], data_line[8], data_line[9]))




