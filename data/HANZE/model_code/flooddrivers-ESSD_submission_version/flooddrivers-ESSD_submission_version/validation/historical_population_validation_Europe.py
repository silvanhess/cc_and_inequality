import geopandas
import numpy as np
import os, sys, inspect
from rasterstats import zonal_stats

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

# Load LAUs
LAU = geopandas.read_file(file('Europe_population_per_LAU'))
NUTS3 = geopandas.read_file(file('NUTS3_regions_from_LAUs'))
countries = np.unique(NUTS3['Country'])

# years
years = list(range(1960, 2020, 10))

# load rasters
pop_2011 = file('Population_100m_baseline')
NUTS_gridpop_2011 = zonal_stats(NUTS3, pop_2011, stats=['sum'], all_touched=True)

Results_all = np.zeros([len(years),80])

for iyy, year in enumerate(years):

    pop_year = file('Population_100m_year')+str(year)+'.tif'
    NUTS_gridpop_year = zonal_stats(NUTS3, pop_year, stats=['sum'], all_touched=True)
    Result_year = np.zeros([len(LAU), 10])

    # zonal statistics
    for idx, nuts in enumerate(NUTS3.itertuples()):

        LAU_NUTS = LAU.loc[LAU['NUTS3']==nuts.NUTS3,]
        if NUTS_gridpop_2011[idx]['sum'] == None:
            continue
        POP_2011_grid = NUTS_gridpop_2011[idx]['sum'] / 100
        POP_year_grid = NUTS_gridpop_year[idx]['sum'] / 100
        POP_2011_table = nuts.SUM_POP_24
        POP_year_table = nuts[iyy+2]
        correction2011 = POP_2011_table / POP_2011_grid
        correctionyear = POP_year_table / POP_year_grid
        LAU_gridpop_2011 = zonal_stats(LAU_NUTS, pop_2011, stats=['sum'], all_touched=True)
        LAU_gridpop_year = zonal_stats(LAU_NUTS, pop_year, stats=['sum'], all_touched=True)
        for idy, lau in enumerate(LAU_NUTS.itertuples()):
            pop_year_c = iyy+11
            if LAU_gridpop_2011[idy]['sum'] == None:
                continue
            POP_2011_grid_lau = LAU_gridpop_2011[idy]['sum'] * correction2011 / 100
            POP_year_grid_lau = LAU_gridpop_year[idy]['sum'] * correctionyear / 100
            if (POP_2011_grid_lau <=0) | (lau.POP_2011 <=0) | (lau[pop_year_c] <=0):
                continue
            else:
                modelled_change = POP_year_grid_lau / POP_2011_grid_lau
                observed_change = lau[pop_year_c] / lau.POP_2011
                modelled_popyear = lau.POP_2011 * modelled_change
                intersect_of_population = (modelled_popyear if modelled_popyear < lau[pop_year_c] else lau[pop_year_c])
                sum_of_population = (modelled_popyear if modelled_popyear > lau[pop_year_c] else lau[pop_year_c])
                HitRate = intersect_of_population / lau[pop_year_c]
                SuccessIndex = intersect_of_population / sum_of_population
                AvgError = abs(modelled_popyear - lau[pop_year_c]) / lau.POP_2011

                Result_year[lau.OBJECTID-1, :] = np.array([lau.OBJECTID, POP_year_grid_lau, POP_2011_grid_lau,
                                                      lau[pop_year_c], lau.POP_2011, modelled_change, observed_change,
                                                      HitRate, SuccessIndex, AvgError])
        print(nuts.NUTS3)

    np.savetxt(file('Europe_validation_results_per_LAU')+str(year)+'.txt', Result_year, fmt='%10.4f', delimiter=',',)
    print(year)

    AbsError = np.absolute(Result_year[:, [4]] * Result_year[:, [5]] - Result_year[:, [3]])
    AbsChange = np.absolute(Result_year[:, [4]] - Result_year[:, [3]])
    AvgError = sum(AbsError) / sum(Result_year[:, [3]])
    AvgRelError = sum(AbsError) / sum(AbsChange)
    Results_all[iyy, 0] = AvgError
    Results_all[iyy, 40] = AvgRelError

    for iyc, c in enumerate(countries):
        ix = LAU['NUTS3'].str.contains(c)
        AbsError_c = np.absolute(Result_year[ix, [4]] * Result_year[ix, [5]] - Result_year[ix, [3]])
        AbsChange_c = np.absolute(Result_year[ix, [4]] - Result_year[ix, [3]])
        AvgError_c = sum(AbsError_c) / sum(Result_year[ix, [3]])
        AvgRelError_c = sum(AbsError_c) / sum(AbsChange_c)
        Results_all[iyy, iyc+1] = AvgError_c
        Results_all[iyy, iyc+41] = AvgRelError_c

np.savetxt(file('Europe_validation_results_aggregate'), Results_all, fmt='%10.4f', delimiter=',',)
