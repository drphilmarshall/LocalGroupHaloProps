# ======================================================================

import numpy as np

# ======================================================================

def draw(x,N):
    
    values = np.zeros(N)
    
    if x[2] == 'Gaussian':
        values += x[0] + x[1]*np.random.randn(N)
    elif x[2] == 'Delta':
        values += x[0]
    else:
        print "ERROR: unrecognised PDF name ",x[2]
        exit()
    
    return values
    
# ======================================================================
