import geopandas
import os, sys, inspect
import pandas as pd
import numpy as np
import rasterio
import re
from rasterio.windows import Window
import matplotlib.pyplot as plt
from matplotlib import colors

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

## this script maps exposure per river and coastal floods from HANZE v1.0
## requires all exposure maps to be computed or downloaded!
Flood_events = pd.read_excel(open(file('Flood_events'), 'rb'), sheet_name='List', index_col='No.')

# Load NUTS regions
NUTS2010 = geopandas.read_file(file('NUTS2010_shapefile'))

for kf, f in enumerate(Flood_events.itertuples()):

    f_year = f.Year
    if f_year<1950:
        f_year=int(np.round(f_year/10)*10)
    elif f_year < 2000:
        f_year = int(np.round(f_year / 5) * 5)
    years = list([1870, f_year, 2020])
    # check type of event
    type = f.Type
    if type=='Coastal':
        hazard_mask_dataset = rasterio.open(file('coastal_flood_mask'))
    elif (type=='River'):
        hazard_mask_dataset = rasterio.open(file('river_flood_mask'))
    else:
        continue

    regions = re.split(";",f.Regions)
    EXT_MIN_Xa = np.zeros([len(regions),1])
    EXT_MAX_Xa = np.zeros([len(regions),1])
    EXT_MIN_Ya = np.zeros([len(regions),1])
    EXT_MAX_Ya = np.zeros([len(regions),1])

    for kr, NUTS_region in enumerate(regions):
        NUTS3_regions_s = NUTS2010.loc[NUTS2010['Code'].str.contains(NUTS_region)]
        # NUTS region codes
        region = NUTS3_regions_s.index.values[0]
        # extent of NUTS region
        EXT_MIN_Xa[kr,0] = NUTS2010['EXT_MIN_X'][region]
        EXT_MAX_Xa[kr,0] = NUTS2010['EXT_MAX_X'][region]
        EXT_MIN_Ya[kr,0] = NUTS2010['EXT_MIN_Y'][region]
        EXT_MAX_Ya[kr,0] = NUTS2010['EXT_MAX_Y'][region]

    # extent of NUTS region
    EXT_MIN_X = EXT_MIN_Xa.min()
    EXT_MAX_X = EXT_MAX_Xa.max()
    EXT_MIN_Y = EXT_MIN_Ya.min()
    EXT_MAX_Y = EXT_MAX_Ya.max()

    start_grid_x = (EXT_MIN_X - 2636000) / 100 - 50
    start_grid_y = 40300 - (EXT_MAX_Y - 1386000) / 100 - 50
    extent_x = (EXT_MAX_X - EXT_MIN_X) / 100 + 100
    extent_y = (EXT_MAX_Y - EXT_MIN_Y) / 100 + 100

    if extent_y > (extent_x * (1 / 1.2)):
        start_grid_x -= int((extent_y * 1.2 - extent_x) * 0.5)
        extent_x = int(extent_y * 1.2)

    if extent_x > (extent_y * 1.2):
        start_grid_y -= int((extent_x / 1.2 - extent_y) * 0.5)
        extent_y = int(extent_x / 1.2)

    fig, axs = plt.subplots(3, 3, constrained_layout=True, figsize=(16, 12))
    fig.suptitle(f.Country + " (" + str(f.Year) + ")", fontsize=16)
    sp = 330
    fmax=1E6
    pmax=150
    # load raster data
    for year in years:
        sp = sp + 1
        name_suffix = str(year) + '.tif'
        clc_dataset = rasterio.open(file('Land_cover_year')+name_suffix)
        pop_dataset = rasterio.open(file('Population_100m_year')+name_suffix)
        fas_dataset = rasterio.open(file('Wealth_year')+name_suffix)
        pop_year = pop_dataset.read(1, window=Window(start_grid_x, start_grid_y, extent_x, extent_y))
        clc_year = clc_dataset.read(1, window=Window(start_grid_x, start_grid_y, extent_x, extent_y))
        fas_year = fas_dataset.read(1, window=Window(start_grid_x, start_grid_y, extent_x, extent_y))
        hazard_mask = hazard_mask_dataset.read(1, window=Window(start_grid_x, start_grid_y, extent_x, extent_y))
        hazard_mask[hazard_mask>1]=0

        if sp == 331:
            if fas_year.max() < 5E6:
                fmax = fas_year.max()/5
            if pop_year.max() < 750:
                pmax = pop_year.max()/5

        # make a color map of fixed colors
        cmap = colors.ListedColormap(['grey','red','purple','yellow','green','blue'])
        bounds=[0,1,3,12,23,41,45]
        norm = colors.BoundaryNorm(bounds, cmap.N)

        # CLC
        plt.subplot(sp)
        im1 = plt.imshow(clc_year, cmap=cmap, norm=norm)
        im1.axes.set_title(str(year), fontsize=11)
        if sp == 331:
            im1.axes.set_ylabel('Land cover/use')
            cbar = plt.colorbar(ax=im1.axes, cmap=cmap, norm=norm, boundaries=bounds, ticks=[0, 1, 3, 12, 23, 41, 45],
                                location='left')
            cbar.ax.set_yticklabels(['','Sea and ocean','Urban','Other artificial','Agriculture','Natural','Inland water'])
        plt.imshow(hazard_mask, cmap='Blues',alpha=0.25)

        # Pop
        plt.subplot(sp+3)
        im2 = plt.imshow(pop_year, cmap='YlOrRd', vmax=pmax)
        if sp == 331:
            im2.axes.set_ylabel('Population density per ha')
            plt.colorbar(ax=im2.axes, location='left')
        plt.imshow(hazard_mask, cmap='Blues',alpha=0.25)

        # FA
        plt.subplot(sp+6)
        im3 = plt.imshow(fas_year, cmap='PuRd', vmax=fmax)
        if sp == 331:
            im3.axes.set_ylabel('Fixed assets per ha (euro)')
            plt.colorbar(ax=im3.axes, location='left')
        plt.imshow(hazard_mask, cmap='Blues',alpha=0.25)

    plt.savefig(file('Event_exposure_fig')+str(f.Index)+'_'+f.Country+'_'+str(f.Year)+'.png', dpi=500)
    plt.close(fig)
