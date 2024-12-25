import numpy as np
import os, sys, inspect
import pandas as pd

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

# load samples of Corine Land Cover and its changes, generated in ArcGIS
CLC_changes = np.genfromtxt(file('CLC_changes_sample_data'), skip_header=True, delimiter=',')
CLC_nochanges = np.genfromtxt(file('CLC_nochanges_sample_data'), skip_header=True, delimiter=',')

CLC_changes2 = np.concatenate((CLC_changes[:,1:3],CLC_changes[:,3:26]),axis=1)
CLC_nochanges2 = np.concatenate((CLC_nochanges[:,[25]],CLC_nochanges[:,[25]],CLC_nochanges[:,2:25]),axis=1)
CLC_data = np.concatenate((CLC_changes2,CLC_nochanges2),axis=0)

# correct urban distances greater than
CLC_data[CLC_data[:,3]<0,[3]] = 2500
CLC_data[CLC_data[:,4]<0,[4]] = 1500
CLC_data[CLC_data[:,5]<0,[5]] = 1500
CLC_data[CLC_data[:,6]<0,[6]] = 1000
CLC_data[CLC_data[:,7]<0,[7]] = 500
CLC_data[CLC_data[:,11]<0,[11]] = 0

ix = np.any(CLC_data < 0, axis=1)
CLC_data = CLC_data[~ix,:]

clc_reclass = np.array([[1, 2, 1],
                        [3, 11, 2],
                        [12, 17, 3],
                        [18, 18, 4],
                        [19, 22, 3],
                        [23, 29, 5],
                        [30, 31, 0],
                        [32, 32, 5],
                        [33, 34, 0],
                        [35, 36, 5],
                        [37, 44, 0],
                        [111, 112, 1],
                        [121, 142, 2],
                        [211, 223, 3],
                        [231, 231, 4],
                        [241, 244, 3],
                        [311, 324, 5],
                        [331, 332, 0],
                        [333, 333, 5],
                        [334, 335, 0],
                        [411, 412, 5],
                        [421, 523, 0],
                        ])

cols = range(0,len(clc_reclass))
for c in cols:
    for d in [0,1]:
        ix = (CLC_data[:,d] >=clc_reclass[c,0]) & (CLC_data[:,d] <=clc_reclass[c,1])
        CLC_data[ix, [d]] = clc_reclass[c,2]

ix = np.any(CLC_data[:,0:2] == 0, axis=1)
CLC_data = CLC_data[~ix,:]

# save relevant data as a Numpy dataset
np.save(file('BN_sample_data'), CLC_data[:,[0,1,9,13,23]])