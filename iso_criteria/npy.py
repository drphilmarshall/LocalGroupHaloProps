import numpy as np
import os
import sys
sys.path.append('/u/ki/yymao/pyscripts')
from helpers.SimulationAnalysis import readHlist

base_path = '/nfs/slac/g/ki/ki20/cosmo/behroozi/Consuelo/'

box = int(sys.argv[1])

#os.mkdir('/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/'+str(box))
halos = readHlist(base_path+str(box)+'/hlists/hlist_1.00000.list')
np.save('/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/'+str(box)+'/'+str(box)+'hlist.npy', halos)
