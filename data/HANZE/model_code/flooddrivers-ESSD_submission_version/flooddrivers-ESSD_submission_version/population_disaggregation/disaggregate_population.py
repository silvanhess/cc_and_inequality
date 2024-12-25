# Description: disaggregate 1 km population grid to a 100 meter resolution

# Import system modules
from rasterio.windows import Window
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import rasterio

import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

def disaggregate_population():

    # define output file
    output_file = file('Population_100m_baseline')

    # Define power function
    def power_func(x, A, b):
        return A * x ** b

    # Read table with data
    df_buildings = pd.read_csv (file('ESM_GEOSTAT_statistics'), delimiter=';')
    df_streets = pd.read_csv (file('ESM_Street_GEOSTAT_statistics'), delimiter=';')
    df_imp = pd.read_csv (file('IMP_GEOSTAT_statistics'), delimiter=';')

    # define data for fitting
    xdata = df_imp['VALUE'].to_numpy()
    data_buildings = df_buildings['MEAN'].to_numpy()
    data_streets = df_streets['MEAN'].to_numpy()
    data_imp = df_imp['MEAN'].to_numpy()

    # fit power function to the data
    popt_buildings, _ = curve_fit(power_func, xdata[1:65], data_buildings[1:65], bounds=(1, [100., 2.]))
    popt_streets, _ = curve_fit(power_func, xdata[1:17], data_streets[1:17], bounds=(1, [200., 2.]))
    popt_imp, _ = curve_fit(power_func, xdata[1:85], data_imp[1:85], bounds=(1, [100., 2.]))

    # load population raster and open other datasets without loading
    geostat_dataset = rasterio.open(file('GEOSTAT_1km_population'))
    geostat_population = geostat_dataset.read(1)
    geostat_grid = geostat_dataset.shape
    columns = list(range(0,geostat_grid[1]))
    rows = list(range(0,geostat_grid[0]))

    clc_dataset = rasterio.open(file('Land_cover_baseline'))
    clc_grid = clc_dataset.shape
    building_dataset = rasterio.open(file('ESM2012_buildings'))
    imp_dataset = rasterio.open(file('ESM2012_streets'))
    street_dataset = rasterio.open(file('Imperviousness2012'))

    # preallocate 100-meter population grid
    pop100 = np.zeros([clc_grid[0],clc_grid[1]], dtype = np.uint16)

    # load population thresholds
    thresholds = np.genfromtxt(file('Population_thresholds'),dtype=float,delimiter=',',skip_header=True)

    # initialize random generator
    rng = np.random.default_rng(12345)

    # disaggregation loop
    for c in columns:
        for r in rows:
            pop1km = geostat_population[r, c]
            if pop1km not in [0,65535]:
                # read raster data only for the particular grid cell
                clc = clc_dataset.read(1, window=Window(c * 10, r * 10, 10, 10))
                clc[clc>44] = 0 # correct for 'not classified' value of 48
                build = building_dataset.read(1, window=Window(c * 10, r * 10, 10, 10))
                imp = imp_dataset.read(1, window=Window(c * 10, r * 10, 10, 10))
                imp[imp>100] = 0 # correct for NoData value of 255
                street = street_dataset.read(1, window=Window(c * 10, r * 10, 10, 10))

                # mask out areas with no CLC class, but some artificial surfaces.
                build = np.where(clc > 0, build, 0)
                imp = np.where(clc > 0, imp, 0)
                street = np.where(clc > 0, street, 0)

                # remove NoData from CLC cells:
                clc_land = clc[clc>0]
                # check if there are any cells. If not, skip and go to next GEOSTAT cell.
                if any(clc_land):
                    # identify CLC class and count their frequency
                    clc_classes, clc_count = np.unique(clc_land, return_counts=True)
                    # find thresholds applicable to CLC classes
                    threshold_local = thresholds[np.isin(thresholds[:,0],clc_classes),:]
                    # sort thresholds by population density (ascending)
                    sorted_thresholds = threshold_local[threshold_local[:,1].argsort(),0]
                    # check if there are any buildings. If not, use imperviousness, and if still nothing, use streets,
                    # if completely no artificial surfaces
                    if build.max()>0:
                        art_surface = build.astype('float')
                        popt_surface = popt_buildings
                    elif imp.max()>0:
                        art_surface = imp.astype('float')
                        popt_surface = popt_imp
                    elif street.max()>0:
                        art_surface = street.astype('float')
                        popt_surface = popt_streets
                    else:
                        art_surface = np.ones([10,10])
                        popt_surface = np.ones([2,1])
                    # preallocate population grid and surplus population
                    local_pop = np.zeros([10,10])
                    local_pop[clc > 0] = pop1km / sum(sum(clc > 0))
                    surplus_pop = 0
                    remaining_area = 100
                    lc = 0
                    # disaggregation loop by CLC patch
                    for l in list(sorted_thresholds):
                        lc = lc + 1
                        # total population assigned to the CLC patch
                        base_population = sum(local_pop[clc == l])
                        # extract artificial surfaces for the CLC patch
                        clc_art_surface = np.where(clc == l, art_surface, 0)
                        # count size of the CLC patch
                        number_cells = clc_count[clc_classes == l]
                        remaining_area = remaining_area - number_cells
                        remaining_clc_classes = sorted_thresholds[lc:]
                        # if there are no artificial surfaces, assign 0 population
                        if sum(sum(clc_art_surface)) == 0:
                            surplus_pop = base_population
                            local_pop[clc == l] = 0
                        else:
                            # maximum population possible in a given CLC patch (rounded to nearest integer up)
                            max_pop = threshold_local[threshold_local[:,0] == l,1]*number_cells/100
                            # if the maximum is lower than the assigned population, adjust population
                            if max_pop < base_population:
                                factor = max_pop / number_cells
                                local_pop[clc == l] = max_pop / number_cells
                                surplus_pop = base_population - max_pop
                            # distribute population according to degree of surface sealing.
                            pop_by_art_surface = clc_art_surface ** popt_surface[1] * popt_surface[0] / 8000
                            factor = sum(local_pop[clc == l])/sum(sum(pop_by_art_surface))
                            local_pop[clc == l] = pop_by_art_surface[clc == l] * factor
                        # compute surplus population
                        if surplus_pop > 0:
                            local_pop = np.where(np.isin(clc, remaining_clc_classes),
                                                 local_pop + surplus_pop / remaining_area,
                                                 local_pop)
                            surplus_pop = 0

                    # if there is still population remaining, distribute proportionally to population already
                    # distributed
                    if sum(sum(local_pop)) < pop1km:
                        adjustment = pop1km / sum(sum(local_pop))
                        local_pop = local_pop * adjustment

                    # Round the population to integer and check if the population matches
                    local_pop_rounded = local_pop.round(decimals=0)
                    if sum(sum(local_pop_rounded)) > pop1km:
                        while sum(sum(local_pop_rounded)) > pop1km:
                            # remove the smallest number of (unrounded) population to make it go to lower integer
                            mod = local_pop % 1 - 0.5
                            diff = mod[mod >= 0].min()
                            local_pop = np.where(mod==diff,local_pop - (diff + 0.000001),local_pop)
                            local_pop_rounded = local_pop.round(decimals=0)
                        pop_diff = sum(sum(local_pop_rounded)) - pop1km
                        # if the resulting population is smaller, add population to the cells with most population
                        if pop_diff < 0:
                            local_pop_values, local_pop_count = np.unique(local_pop, return_counts=True)
                            d = len(local_pop_values)
                            while pop_diff < 0:
                                d = d - 1
                                # check if there are more cells with equal population than the missing population
                                if local_pop_count[d] > -pop_diff:
                                    value_indices = np.argwhere(local_pop == local_pop_values[d])
                                    # if there are, see if they can be distinguished by imperviousness
                                    imp_values = list(np.unique(imp[local_pop == local_pop_values[d]]))
                                    street_values = list(np.unique(street[local_pop == local_pop_values[d]]))
                                    if len(imp_values) >= local_pop_count[d]:
                                        imp_values.sort(reverse=True)
                                        for i in imp_values[:int(-pop_diff)]:
                                            local_pop = np.where((local_pop == local_pop_values[d])&(imp==i),
                                                                 local_pop+1,local_pop)
                                    # if not, see if they can be distinguished by streets
                                    elif len(street_values) >= local_pop_count[d]:
                                        street_values.sort(reverse=True)
                                        for s in street_values[:int(-pop_diff)]:
                                            local_pop = np.where((local_pop == local_pop_values[d])&(street==s),
                                                                 local_pop+1,local_pop)
                                    # if not, randomly assign values
                                    else:
                                        locs = list(rng.integers(low=0, high=local_pop_count[d], size=int(-pop_diff)))
                                        for l in locs:
                                            local_pop[value_indices[l,0],value_indices[l,1]] += 1
                                else:
                                    local_pop = np.where(local_pop == local_pop_values[d], local_pop + 1, local_pop)

                                local_pop_rounded = local_pop.round(decimals=0)
                                pop_diff = sum(sum(local_pop_rounded)) - pop1km

                    elif sum(sum(local_pop_rounded)) < pop1km:
                        while sum(sum(local_pop_rounded)) < pop1km:
                            # remove the smallest number of (unrounded) population to make it go to lower integer
                            mod = local_pop % 1 - 0.5
                            diff = mod[mod <= 0].max()
                            local_pop = np.where(mod==diff, local_pop + (0.500001 - diff), local_pop)
                            local_pop_rounded = local_pop.round(decimals=0)
                        pop_diff = sum(sum(local_pop_rounded)) - pop1km
                        # if the resulting population is smaller, add population to the cells with most population
                        if pop_diff > 0:
                            local_pop_values = np.unique(local_pop)
                            d = len(local_pop_values)
                            while pop_diff > 0:
                                d = d - 1
                                local_pop_count = sum(sum(local_pop == local_pop_values[d]))
                                # check if there are more cells with equal population than the missing population
                                if local_pop_count > pop_diff:
                                    value_indices = np.argwhere(local_pop == local_pop_values[d])
                                    # if there are, see if they can be distinguished by imperviousness
                                    imp_values = list(np.unique(imp[local_pop == local_pop_values[d]]))
                                    street_values = list(np.unique(street[local_pop == local_pop_values[d]]))
                                    if len(imp_values) >= local_pop_count:
                                        for i in imp_values[:int(pop_diff)]:
                                            local_pop = np.where((local_pop == local_pop_values[d])&(imp==i),
                                                                 local_pop-1,local_pop)
                                    # if not, see if they can be distinguished by streets
                                    elif len(street_values) >= local_pop_count:
                                        for s in street_values[:int(pop_diff)]:
                                            local_pop = np.where((local_pop == local_pop_values[d])&(street==s),
                                                                 local_pop-1,local_pop)
                                    # if not, randomly assign values
                                    else:
                                        while pop_diff > 0:
                                            loc = rng.integers(low=0, high=local_pop_count, size=int(1))
                                            if local_pop[value_indices[loc,0],value_indices[loc,1]] >= 1:
                                                local_pop[value_indices[loc,0],value_indices[loc,1]] -= 1
                                                pop_diff = pop_diff - 1
                                else:
                                    local_pop = np.where(local_pop == local_pop_values[d], local_pop - 1, local_pop)

                                local_pop_rounded = local_pop.round(decimals=0)
                                pop_diff = sum(sum(local_pop_rounded)) - pop1km

                    # check to be sure that the disaggregation was correct
                    assert sum(sum(local_pop_rounded)) == pop1km

                    pop100[r*10:(r+1)*10,c*10:(c+1)*10] = local_pop_rounded
                    print(str(c)+", "+str(r))

    # save population raster
    with rasterio.Env():

        profile = clc_dataset.profile
        profile.update(dtype=rasterio.uint16,count=1,compress='lzw',nodata = 65535)
        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(pop100, 1)
        print('Output file saved')
