import numpy as np
import pandas as pd
import os, sys, inspect
import geopandas as gp
import re

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

## Aggregate HANZE flood events and compute sum of events per region

# load HANZE data
Events_hanze_v1 = pd.read_excel(open(file('Flood_events_v1'), 'rb'), sheet_name='List', index_col='No.')
Events_hanze_v2 = pd.read_excel(open(file('Flood_events_v2'), 'rb'), sheet_name='List', index_col='ID')
Events_hanze_v2_B_list = pd.read_excel(open(file('Flood_events_v2_B_list'), 'rb'), sheet_name='List', index_col='No.')

# NUTS 2010 with difference from V1

NUTS2010 = gp.read_file(file('NUTS2010_shapefile'))
NUTS_code = NUTS2010['Code'].values

Regions_v1 = Events_hanze_v1['Regions'].values
Regions_v2 = Events_hanze_v2['Regions'].values
Regions_v2_B_list = Events_hanze_v2_B_list['Regions'].values

Reg_count = np.zeros([len(NUTS2010),4],dtype='object')
Reg_count[:,0] = NUTS_code

for R in Regions_v1:
    regions = re.split(";", R)
    for r in regions:
        ix = NUTS_code==r
        Reg_count[ix,1] +=1

for R in Regions_v2:
    regions = re.split(";", R)
    for r in regions:
        ix = NUTS_code==r
        Reg_count[ix,2] +=1
        if sum(ix)==0:
            print(r)

for R in Regions_v2_B_list:
    regions = re.split(";", R)
    for r in regions:
        ix = NUTS_code==r
        if sum(ix)==0:
            print(r)

Reg_count[:,3] = Reg_count[:,2] - Reg_count[:,1]

cols = ['NUTS3','Floods_v1','Floods_V2','Diff']

loss_events_df = pd.DataFrame(data=Reg_count,index=NUTS2010.index, columns=cols)
loss_events_df.to_csv(file('Flood_events_summary_diff'), sep=',')

# NUTS 2021

NUTS2021 = gp.read_file(file('NUTS2021_shapefile'))
NUTS_code2021 = NUTS2021['NUTS'].values

Regions_v2_2021 = Events_hanze_v2['Region2021'].values
Regions_v2_2021_B_list = Events_hanze_v2['Regions2021'].values

Reg_count2021 = np.zeros([len(NUTS2021),2],dtype='object')
Reg_count2021[:,0] = NUTS_code2021

for R in Regions_v2_2021:
    regions = re.split(";", R)
    for r in regions:
        ix = NUTS_code2021==r
        Reg_count2021[ix,1] +=1
        if sum(ix)==0:
            print(r)

for R in Regions_v2_2021_B_list:
    regions = re.split(";", R)
    for r in regions:
        ix = NUTS_code2021==r
        if sum(ix)==0:
            print(r)

cols = ['NUTS3','Floods_V2']

loss_events_df2021 = pd.DataFrame(data=Reg_count2021,index=NUTS2021.index, columns=cols)
loss_events_df2021.to_csv(file('Flood_events_summary_v2'), sep=',')

## Wildfires

main_path = 'C:/HANZE2_temp/EFAS_catchments/'
valid_path = main_path + 'Validation/'

Events_wildfire = pd.read_excel(open(valid_path + 'fire-report.xlsx', 'rb'), sheet_name='wildfires', index_col='ID')
Regions_fire = Events_wildfire['Regions affected'].values

Reg_count_fire = np.zeros([len(NUTS2021),2],dtype='object')
Reg_count_fire[:,0] = NUTS_code2021

for R in Regions_fire:
    regions = re.split(";", R)
    for r in regions:
        ix = NUTS_code2021==r
        Reg_count_fire[ix,1] +=1
        if sum(ix)==0:
            print(r)

cols = ['NUTS3','Fires']

loss_events_df2021 = pd.DataFrame(data=Reg_count_fire,index=NUTS2021.index, columns=cols)
loss_events_df2021.to_csv(valid_path + 'Wildfire_Regions.csv', sep=',')

## Windstorms

Events_windstorm = pd.read_excel(open(valid_path + 'windstorm-report.xlsx', 'rb'), sheet_name='windstorms', index_col='ID')
Regions_wind2 = Events_windstorm['Regions affected (NUTS2)'].values
Regions_wind3 = Events_windstorm['Regions affected (NUTS3)'].values

Reg_count_wind = np.zeros([len(NUTS2021), 2], dtype='object')
Reg_count_wind[:, 0] = NUTS_code2021

for k,R in enumerate(Regions_wind2):
    regions2 = re.split(";", R)
    regions3s = Regions_wind3[k]
    if not isinstance(regions3s,str):
        regions3s = 'XX999'
    regions3 = np.array(re.split(";", regions3s))
    for r in regions2:
        if r in regions3s:
            iy = [r in i for i in regions3]
            regions3_within2 = regions3[np.array(iy)]
            for r3 in regions3_within2:
                ix = NUTS_code2021 == r3
                Reg_count_wind[ix, 1] += 1
                if sum(ix) == 0:
                    print(r)
        else:
            ix = [r in i for i in NUTS_code2021]
            Reg_count_wind[np.array(ix), 1] +=1
            if sum(np.array(ix))==0:
                print(r)

cols = ['NUTS3','Storms']

loss_events_df2021 = pd.DataFrame(data=Reg_count_wind, index=NUTS2021.index, columns=cols)
loss_events_df2021.to_csv(valid_path + 'Windstorms_Regions.csv', sep=',')