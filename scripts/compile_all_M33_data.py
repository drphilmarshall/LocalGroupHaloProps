import numpy as np
import sys

sys.path.append('/u/ki/yymao/pyscripts')
from helpers.SimulationAnalysis import readHlist

header = '# idMW idM31 idM33 Mvir x y z vx vy vz\n'

path = '/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Consuelo/LG_project/'
M33_prefix = 'MW_halos_'
M33_suffix = '_z0p000_vir_bgc2.txt'

M33_file = open('MW_all_z0p000_data.txt','a+')
M33_file.write(header)
for i in range(4001,4051):
    M33_file_name = path+M33_prefix+str(i)+M33_suffix
    M33_data = readHlist(M33_file_name)
    print M33_file_name
    for data_line in M33_data:
       # print (data_line)
        M33_file.write("%d %d %d %f %f %f %f %f %f %f\n"%(data_line[0], data_line[1], data_line[2], data_line[3], data_line[4], data_line[5], data_line[6], data_line[7], data_line[8], data_line[9]))

M33_file.close()
