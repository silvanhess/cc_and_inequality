# Description: use pure population-land cover (GEOSTAT - Corine Land Cover) cells and compute 70th percentile of
# population density per 1 km^2.

# Import system modules
import pandas as pd
import numpy as np
import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

# Import CSV file with population in pure CLC 1 km grid cells
df = pd.read_csv (file('Population_1km_CLC_pure_cells'), delimiter=';')

# Aggregate CLC classes 122-142
df['U2018_CLC2'] = df['U2018_CLC2'].replace(list(range(4,12)),99)

# Unique CLC values occurring in the dataset
clc = pd.unique(df['U2018_CLC2'])

# Calculate 70th percentile
thresholds = pd.DataFrame(np.zeros((44, 1)),columns=['thresholds'],index=list(range(1,45)))  #np.full((44, 1), 0)
uninhabitable_clc = np.concatenate(([30,31], np.array(range(33,45))), axis=0)
for c in clc:
    if c not in uninhabitable_clc:
        df_c = df.loc[df['U2018_CLC2']==c,'RASTERVALU']
        threshold = np.percentile(df_c,70,interpolation='linear')
        if c==99:
            thresholds.loc[range(4,12),'thresholds'] = round(threshold)
        else:
            thresholds.loc[c,'thresholds'] = round(threshold)

# set dummy threshold to yninhabitable CLCs to avoid divide-by-zero problem
thresholds.loc[uninhabitable_clc,'thresholds'] = 0.0001

# Export table
thresholds.to_csv(file('Population_thresholds'))