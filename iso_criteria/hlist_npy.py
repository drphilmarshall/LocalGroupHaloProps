import os
import numpy as np

#boxes = [4001, 4002, 4003, 4004, 4020, 4026, 4027, 4028, 4029, 4030, 4032, 4033, 4034, 4035, 4036, 4039]
boxes = [4001]
base_path = '/nfs/slac/g/ki/ki20/cosmo/behroozi/Consuelo/'

for box in boxes:
    os.system('bsub -q long python npy.py '+str(box))
