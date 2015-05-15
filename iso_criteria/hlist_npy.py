import os
import numpy as np

#boxes = [4026, 4027, 4028, 4029, 4030, 4032, 4033, 4034, 4035, 4036]
boxes = [4035, 4036]
base_path = '/nfs/slac/g/ki/ki20/cosmo/behroozi/Consuelo/'

for box in boxes:
    os.system('bsub -q kipac-ibq python npy.py '+str(box))
