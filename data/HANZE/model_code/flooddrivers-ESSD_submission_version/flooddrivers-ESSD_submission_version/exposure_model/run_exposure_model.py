from exposure_model import main_model
import concurrent.futures

### RUN EXPOSURE MODEL
## remember to define path to data in get_file.py
## define parameters
# timespan of the simulation (must be an iterable and must match the input NUTS3 database)
years = list(range(2020, 2021)) #list(range(1870, 1950, 10)) + list(range(1950, 2000, 5)) +list(range(2000, 2021))

# regions for which run the model; it will try find a match in the string (doesn't have to be the full string) among
# NUTS3 codes (e.g. 'DE600' for Hamburg). If it doesn't find a match (e.g. "ALL" is written), all regions will be used
regions_option = range(0,1422)

# optional: define a raster (from get_file script) with hazard zones, for which to carry out the calculation
# it is required for uncertainty analysis (provide non-string to skip)!
# insert 'coastal_flood_mask' for 100-year coastal flood hazard zone
# insert 'river_flood_mask' for 100-year riverine flood hazard zone
# if option not needed, provide an input other than string, e.g. "0"
hazard_mask = 0

# optional: do an uncertainty analysis. Requires defining number of iterations (e.g. 100).
# if not desired, insert "0"
uncertainty_option = 0

# optional: compute flood damages using JRC damage functions and water depths (flood footprints only a.t.m.)
# doesn't support yet uncertainty analysis
compute_damage = 0

## run the model (single thread)
# main_model(years, regions_option, hazard_mask, uncertainty_option, compute_damage)

## Note: raster maps per variable are saved, unless the hazard or hazard/uncertainty option is selected.
## In that case, text files with uncertainty bounds per variable and region

## run the model (for parallel threads)
# for generation of exposure maps, parallelize by year.
if hazard_mask == 0:
    cpu_factor = 1
    iterator = years
# for generation of exposure estimates per hazard zone, parallelize by NUTS3 region.
else:
    cpu_factor = 8
    iterator = regions_option

if len(iterator)>=(4 * cpu_factor):
    threads = 2 * cpu_factor
else:
    threads = len(iterator)

def parallel_run(i):
    try:
        if hazard_mask == 0:
            main_model(list([i]), regions_option, 0, 0, 0)
        else:
            main_model(years, i, hazard_mask, uncertainty_option, compute_damage)
    except:
        print('error with %d' % i)

def main():
    executor = concurrent.futures.ProcessPoolExecutor(threads)
    [executor.submit(parallel_run, i) for i in iterator]

if __name__ == '__main__':
    main()