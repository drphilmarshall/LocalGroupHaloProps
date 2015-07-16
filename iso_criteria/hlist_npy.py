import os
import numpy as np

#boxes = [4002, 4003, 4004, 4020, 4026, 4027, 4028, 4029, 4030, 4032, 4033, 4034, 4035, 4036, 4039, 4040]
boxes = [4035, 4036, 4039, 4040]

for box in boxes:
    os.system('bsub -q kipac-ibq python npy.py '+str(box))
