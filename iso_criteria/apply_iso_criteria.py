import iso_criteria
import numpy as np
import sys
sys.path.append('/u/ki/yymao/pyscripts')
from fast3tree import fast3tree

vmax = int(sys.argv[1])
base_path = sys.argv[2]
box_num = int(sys.argv[3])

hlist_path = base_path+str(box_num)+'/'+str(box_num)+'hlist.npy'
if vmax:
    vmax_cut=vmax
else:
    vmax_cut=None

halos = np.load(hlist_path)
halos = halos[halos['vmax']>80]
host_flag = np.where(halos['upid'] == -1)[0]

pair_dist = 1.0/0.7 #units of Mpc/h
M33_dist = 0.4/0.7 #units of Mpc/h
iso_dist = 2.0/0.7 #units of Mpc/h
box_size = 420 #units of Mpc/h, consuelo

pairs = iso_criteria.criteria.find_pairs(halos, host_flag, pair_dist)
MW_M31_larger, M31_M31_larger, M33_M31_larger, M33_M33_larger, M31_M33_larger, MW_M33_larger = iso_criteria.criteria.isolated(pairs, halos, host_flag, iso_dist, M33_dist, box_size, vmax_cut, periodic=True)

dat_M31_larger = iso_criteria.criteria.get_data(halos, host_flag, MW_M31_larger, M31_M31_larger, M33_M31_larger)
dat_M33_larger = iso_criteria.criteria.get_data(halos, host_flag, MW_M33_larger, M31_M33_larger, M33_M33_larger)
trip_dat_M31_larger = iso_criteria.criteria.get_trip_data(dat_M31_larger)
trip_dat_M33_larger = iso_criteria.criteria.get_trip_data(dat_M33_larger)

save_path = base_path+str(box_num)+'/'
if vmax:
    np.save(save_path+'dat_M31_larger_vmax>'+str(vmax)+'.npy', dat_M31_larger)
    np.save(save_path+'dat_M33_larger_vmax>'+str(vmax)+'.npy', dat_M33_larger)
    np.save(save_path+'trip_dat_M31_larger_vmax>'+str(vmax)+'.npy', trip_dat_M31_larger)
    np.save(save_path+'trip_dat_M33_larger_vmax>'+str(vmax)+'.npy', trip_dat_M33_larger)
else:
    np.save(save_path+'dat_M31_larger_vres80.npy', dat_M31_larger)
    np.save(save_path+'dat_M33_larger_vres80.npy', dat_M33_larger)
    np.save(save_path+'trip_dat_M31_larger_vres80.npy', trip_dat_M31_larger)
    np.save(save_path+'trip_dat_M33_larger_vres80.npy', trip_dat_M33_larger)
