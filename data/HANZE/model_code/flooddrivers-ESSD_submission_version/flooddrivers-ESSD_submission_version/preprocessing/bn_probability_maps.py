from pomegranate import *
from rasterio.windows import Window
import geopandas
import numpy as np
import os, sys, inspect
import pandas as pd
import rasterio

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file
from auxiliary_functions import load_dataset_by_NUTS, NUTS_mask, write_empty_raster

# prepare BN
CLC_data = np.load(file('BN_sample_data'))

VLAU_bins = pd.qcut(CLC_data[:, 2], 9)
VLAU = DiscreteDistribution.from_samples(VLAU_bins.codes)
VLAU_data = np.expand_dims(VLAU_bins.codes, axis=1)
ylHr0_whe_bins = pd.qcut(CLC_data[:, 4], 5)
ylHr0_whe = DiscreteDistribution.from_samples(ylHr0_whe_bins.codes)
ylHr0_whe_data = np.expand_dims(ylHr0_whe_bins.codes, axis=1)
gras200a_yld_bins = pd.qcut(CLC_data[:, 3], 10)
gras200a_yld = DiscreteDistribution.from_samples(gras200a_yld_bins.codes)
gras200a_yld_data = np.expand_dims(gras200a_yld_bins.codes, axis=1)

BN_data = np.concatenate((VLAU_data, ylHr0_whe_data, gras200a_yld_data, CLC_data[:, [0]]), axis=1)
CLC_before = ConditionalProbabilityTable.from_samples(BN_data, [VLAU, ylHr0_whe, gras200a_yld])

BN_data = np.concatenate((CLC_data[:, [0]], VLAU_data, ylHr0_whe_data, gras200a_yld_data, CLC_data[:, [1]]), axis=1)
cpt = ConditionalProbabilityTable.from_samples(BN_data, [CLC_before, VLAU, ylHr0_whe, gras200a_yld])

BN_data = np.concatenate((VLAU_data, ylHr0_whe_data, gras200a_yld_data, CLC_data[:, [1]]), axis=1)
CLC_after2 = ConditionalProbabilityTable.from_samples(BN_data, [VLAU, ylHr0_whe, gras200a_yld])

BN_data = np.concatenate((CLC_data[:, [1]], VLAU_data, ylHr0_whe_data, gras200a_yld_data, CLC_data[:, [0]]), axis=1)
cpt2 = ConditionalProbabilityTable.from_samples(BN_data, [CLC_after2, VLAU, ylHr0_whe, gras200a_yld])

# define raster datasets
nuts_dataset = rasterio.open(file('NUTS2010_raster'))
pop_dataset = rasterio.open(file('Population_100m_baseline'))
clc_dataset = rasterio.open(file('Land_cover_baseline'))
virtual_laus_pop_den = rasterio.open(file('Population_density_vLAU'))
wheat_dataset = rasterio.open(file('GAEZ_wheat'))
grass_dataset = rasterio.open(file('GAEZ_grass'))
bn_maps = file('BN_map')

# Load NUTS regions
NUTS2010 = geopandas.read_file(file('NUTS2010_shapefile'))
NUTS_code = NUTS2010['Code']
NUTS3_regions = pd.DataFrame(NUTS_code)

params = ['to_urban','to_crop','to_past','from_crop','from_past','to_crop_p','to_past_p','from_crop_p','from_past_p']

# prepare empty rasters
dims = np.zeros([clc_dataset.height,clc_dataset.width], dtype = np.uint16)
for param in params:
    write_empty_raster(param, pop_dataset.profile, bn_maps, np.uint16, dims)

bn_states = range(0,11250)

# regional loop
for nuts3 in NUTS3_regions.itertuples():
    region = nuts3.Index
    region_code = NUTS_code[region]
    # find grid cells specific for the NUTS region
    region_mask, location = NUTS_mask(NUTS2010, region, nuts_dataset)
    # load other raster data
    clc_baseline = load_dataset_by_NUTS(clc_dataset, location, region_mask, 0)
    local_pop_density = load_dataset_by_NUTS(virtual_laus_pop_den, location, region_mask, 0)
    wheat_yield = load_dataset_by_NUTS(wheat_dataset, location, region_mask, 0)
    grass_yield = load_dataset_by_NUTS(grass_dataset, location, region_mask, 0)

    clc_baseline_1d = clc_baseline[region_mask]
    local_pop_density_1d = local_pop_density[region_mask]
    wheat_yield_1d = wheat_yield[region_mask]
    grass_yield_1d = grass_yield[region_mask]

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
                            [37, 50, 0]
                            ])

    cols = range(0, len(clc_reclass))
    for c in cols:
        for d in [0, 1]:
            ix = (clc_baseline_1d >= clc_reclass[c, 0]) & (clc_baseline_1d <= clc_reclass[c, 1])
            clc_baseline_1d[ix] = clc_reclass[c, 2]

    X = range(0, clc_baseline_1d.size)
    s_pop_v = pd.cut(local_pop_density_1d, VLAU_bins.categories).codes
    s_whe_v = pd.cut(wheat_yield_1d, ylHr0_whe_bins.categories).codes
    s_gra_v = pd.cut(grass_yield_1d, gras200a_yld_bins.categories).codes

    classes_from = [[2, 3, 4, 5], [1, 2, 4, 5], [1, 2, 3, 5], [3], [4], [1, 2, 4, 5], [1, 2, 3, 5], [3], [4]]
    classes_to = [[1], [3], [4], [3], [4], [3], [4], [3], [4]]

    for k, p in enumerate(params):
        bn_probability_1d = np.zeros([clc_baseline_1d.size])
        for b in bn_states:
            if k <= 4:
                cond = cpt.parameters[0][b]
            else:
                cond = cpt2.parameters[0][b]
            if (cond[0] in classes_from[k]) & (cond[4] in classes_to[k]):
                ix = (clc_baseline_1d == cond[0]) & (s_pop_v == cond[1]) & (s_whe_v == cond[2]) & (s_gra_v == cond[3])
                if (k == 3) | (k == 4) | (k >= 7):
                    bn_probability_1d[ix] = 1 - cond[5]
                else:
                    bn_probability_1d[ix] = cond[5]

        bn_dataset_1 = rasterio.open(bn_maps+p+'.tif')
        read_dataset_bn_1 = bn_dataset_1.read(1, window=Window(location[0], location[1], location[2], location[3]))
        read_dataset_bn_1[region_mask] = bn_probability_1d * 10000
        profile = bn_dataset_1.profile
        bn_dataset_1.close()
        with rasterio.open(bn_maps+p+'.tif', 'r+', **profile) as dst:
            dst.write(read_dataset_bn_1, window=Window(location[0], location[1], location[2], location[3]), indexes=1)
        bn_dataset_1.close()
        print(region_code)