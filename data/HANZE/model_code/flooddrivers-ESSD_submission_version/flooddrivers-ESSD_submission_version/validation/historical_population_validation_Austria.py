import geopandas
import numpy as np
import os, sys, inspect
from rasterstats import zonal_stats

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

# Load LAUs
LAU = geopandas.read_file(file('Austria_population_per_LAU'))

# years
years = list(range(1870, 1950, 10)) + list(range(1950, 2025, 5))

# compute baseline population per LAU from the 100 m grid
LAU_gridpop_2011 = zonal_stats(LAU, file('Population_100m_baseline'), stats=['sum'], all_touched=True)

Results_all = np.zeros([len(years),4])

for iyy, year in enumerate(years):

    pop_year = file('Population_100m_year')+str(year)+'.tif'
    LAU_gridpop_year = zonal_stats(LAU, pop_year, stats=['sum'], all_touched=True)
    Result_year = np.zeros([len(LAU),10])

    # zonal statistics
    for idy, lau in enumerate(LAU.itertuples()):

        POP_2011_grid_lau = LAU_gridpop_2011[idy]['sum']
        POP_year_grid_lau = LAU_gridpop_year[idy]['sum']
        POP_2011_table = lau.Austria_22
        POP_year_table = lau[iyy+4]
        try:
            if (POP_2011_grid_lau <=0) | (POP_2011_table <=0) | (POP_year_table <=0):
                continue
            else:
                modelled_change = POP_year_grid_lau / POP_2011_grid_lau
                observed_change = POP_year_table / POP_2011_table
                modelled_pop_year = POP_2011_table * modelled_change
                intersect_of_population = (modelled_pop_year if modelled_pop_year < POP_year_table else POP_year_table)
                sum_of_population = (modelled_pop_year if modelled_pop_year > POP_year_table else POP_year_table)
                HitRate = intersect_of_population / POP_year_table
                SuccessIndex = intersect_of_population / sum_of_population
                AvgError = abs(modelled_pop_year - POP_year_table) / POP_year_table
                Result_year[lau.FID, :] = np.array([lau.FID, POP_year_grid_lau, POP_2011_grid_lau, POP_year_table,
                                                    POP_2011_table, modelled_change, observed_change, HitRate,
                                                    SuccessIndex, AvgError])
        except:
            a=1

    np.savetxt(file('Austria_validation_results_per_LAU')+str(year)+'.txt', Result_year, fmt='%10.4f', delimiter=',',)
    print(year)

    AbsError = np.absolute( Result_year[:, [4]]*Result_year[:, [5]] - Result_year[:, [3]] )
    AbsChange = np.absolute(Result_year[:, [4]] - Result_year[:, [3]])
    AvgError = sum(AbsError) / sum(Result_year[:, [3]])
    AvgRelError = sum(AbsError) / sum(AbsChange)
    Results_all[iyy, 0] = AvgError
    Results_all[iyy, 1] = AvgRelError

np.savetxt(file('Austria_validation_results_aggregate'), Results_all, fmt='%10.4f', delimiter=',',)

