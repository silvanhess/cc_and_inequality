import numpy as np
import pandas as pd
import os, sys, inspect
import geopandas as gp
import re

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

## Convert HANZE database into a shapefile

# load HANZE
Events_hanze_v2 = pd.read_excel(open(file('Flood_events_v2'), 'rb'), sheet_name='List')
Events_hanze_v2['geometry'] = np.zeros(len(Events_hanze_v2))
Codes = Events_hanze_v2['ID'].values

Events_hanze_v2['Start_D'] = Events_hanze_v2['Start'].dt.day
Events_hanze_v2['Start_M'] = Events_hanze_v2['Start'].dt.month
Events_hanze_v2['Start_Y'] = Events_hanze_v2['Start'].dt.year
Events_hanze_v2['End_D'] = Events_hanze_v2['End'].dt.day
Events_hanze_v2['End_M'] = Events_hanze_v2['End'].dt.month
Events_hanze_v2['End_Y'] = Events_hanze_v2['End'].dt.year
Events_hanze_v2 = Events_hanze_v2.drop(columns=['Start','End'])

# load NUTS regions
NUTS2010 = gp.read_file(file('NUTS2010_shapefile_s'))
NUTS_code2010 = NUTS2010['Code'].values
NUTS2021 = gp.read_file(file('NUTS2021_shapefile_s'))
NUTS_code2021 = NUTS2021['NUTS'].values

# create empty GDF
Flood_events2010 = gp.GeoDataFrame(columns=Events_hanze_v2.columns, crs="EPSG:3035", geometry='geometry')
Flood_events2021 = gp.GeoDataFrame(columns=Events_hanze_v2.columns, crs="EPSG:3035", geometry='geometry')

# collect all regions from GIS file that belong to each event, dissolve and combine with HANZE data
for C in Codes:
    print(str(C))
    HANZE_data_r = Events_hanze_v2.loc[Events_hanze_v2['ID'] == C,]

    # v2010
    R = HANZE_data_r['Region2010'].values
    regions = np.array(re.split(";", R[0]))
    ix = np.isin(NUTS_code2010,regions)
    NUTS_r = NUTS2010.loc[ix]
    NUTS_r_d = NUTS_r.dissolve()
    NUTS_r_d_g = gp.GeoDataFrame(HANZE_data_r, geometry=NUTS_r_d['geometry'].values)
    Flood_events2010 = pd.concat([Flood_events2010, NUTS_r_d_g], axis=0)

    # v2021
    R = HANZE_data_r['Region2021'].values
    regions = np.array(re.split(";", R[0]))
    ix = np.isin(NUTS_code2021,regions)
    NUTS_r = NUTS2021.loc[ix]
    NUTS_r_d = NUTS_r.dissolve()
    NUTS_r_d_g = gp.GeoDataFrame(HANZE_data_r, geometry=NUTS_r_d['geometry'].values)
    Flood_events2021 = pd.concat([Flood_events2021, NUTS_r_d_g], axis=0)

# Export to shapefile
Flood_events2010.to_file(file('HANZE_SHP_2010'))
Flood_events2021.to_file(file('HANZE_SHP_2021'))