from scipy.stats import rankdata
import numpy as np
import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from auxiliary_functions import copula_inference_normal, copula_inference_Frank

# Subregional (vLAU) population growth
def model_suregional_population(population_baseline, vLAUs, local_pop_density, urban_mask, copula_rural, copula_city,
                                copula_samples, urban_distance_6, copula_sample, uncertainty_mode):

    vLAU_PG = np.zeros([vLAUs.shape[0], 8])
    vLAU_PG[:, 0] = vLAUs
    for idx, vLAU in enumerate(vLAUs):
        PD = vLAU / 100  # population density per km^2 (original data: per 100 km^2)
        if PD < 1500:
            Pop_growth = copula_inference_normal(copula_rural, copula_samples, PD)  # 50-year growth rate
        else:
            KD = np.mean(urban_distance_6[local_pop_density == vLAU])
            Pop_growth = copula_inference_Frank(copula_city, copula_samples, KD)  # 50-year growth rate
        # annual growth rate per subregional unit (vLAU)
        if uncertainty_mode == 1:
            Pop_growth_q = Pop_growth[copula_sample]
        else:
            Pop_growth_q = np.median(Pop_growth)
        Mean_pop_growth = 1 / (Pop_growth_q ** (1 / 50))  # (Pop_growth[idr] ** (1/50))
        if Mean_pop_growth < 0.97743:
            vLAU_PG[idx, 1] = 0.97743
        elif Mean_pop_growth > 1.016464:
            vLAU_PG[idx, 1] = 1.016464
        else:
            vLAU_PG[idx, 1] = Mean_pop_growth
        urban_mask_lau = (local_pop_density == vLAU) & urban_mask
        rural_mask_lau = (local_pop_density == vLAU) & (~urban_mask)
        vLAU_PG[idx, 2] = sum(population_baseline[urban_mask_lau])
        vLAU_PG[idx, 3] = sum(population_baseline[rural_mask_lau])

    return vLAU_PG

# Urban population and urban fabric (CLC111-112)
def model_urban_population(clc_year, population_year, year, region_code, population_baseline, vLAUs,
                           vLAU_PG, local_pop_density, NUTS3_pph, pph_baseline, NUTS3_urban_fraction, NUTS3_population,
                           combined_urban_distance, urban_pop_correction, rural_pop_correction, urban_mask,
                           airports, reservoirs, urban_remove_threshold, urban_add_threshold, landuse_bn_to_urban):


    pph_year = NUTS3_pph[year][region_code]
    urban_pop_year_NUTS = NUTS3_urban_fraction[year][region_code] / 100 * NUTS3_population[year][
        region_code] * 1000 * urban_pop_correction
    rural_pop_year_NUTS = (1 - NUTS3_urban_fraction[year][region_code] / 100) * NUTS3_population[year][
        region_code] * 1000 * rural_pop_correction
    # compute population that should be per vLAU
    if year <= 2011:
        vLAU_PG[:, 4:6] = (vLAU_PG[:, [1]] ** (2011 - year) * vLAU_PG[:, 2:4])
    elif year > 2011:
        vLAU_PG[:, 4:6] = ((1 / vLAU_PG[:, [1]]) ** (year - 2011) * vLAU_PG[:, 2:4])
    vLAU_PG[:, [6]] = vLAU_PG[:, [4]] * (urban_pop_year_NUTS / sum(vLAU_PG[:, [4]]))  # desired urban pop
    vLAU_PG[:, [7]] = vLAU_PG[:, [5]] * (rural_pop_year_NUTS / sum(vLAU_PG[:, [5]]))  # desired rural pop
    # redistribute population within each vLAU
    for vLAU in vLAUs:
        idx = vLAU_PG[:, 0] == vLAU
        if vLAU_PG[idx, [6]] > 1.5:  # check if there is more than one person in urban areas
            urban_mask_lau = (local_pop_density == vLAU) & urban_mask #& pop_baseline_mask
            # urban pop step 1 - modify population according to change in household size
            clc_year_1d = clc_year[urban_mask_lau].flatten()
            population_year[urban_mask_lau] = population_baseline[urban_mask_lau] * (pph_year / pph_baseline)
            population_year_1d = population_year[urban_mask_lau].flatten()
            # step 3 - calculate surplus population
            urban_pop_year = np.rint(vLAU_PG[idx, [6]])
            surplus_urban_pop = np.rint(vLAU_PG[idx, [2]] * pph_year / pph_baseline) - urban_pop_year
            # step 4 - remove surplus population
            if surplus_urban_pop > 0:
                if year <= 2011:
                    udr_unique, udr_indices = np.unique(combined_urban_distance[urban_mask_lau], return_inverse=True)
                    i = udr_unique.size
                    if i > 1:
                        max_distance = np.log((udr_unique[i - 1] - udr_unique[0]))
                    else:
                        max_distance = 1
                    while surplus_urban_pop > 0:
                        i = i - 1
                        loc = udr_indices == i
                        if (udr_unique[i] - udr_unique[0]) <= 1:
                            deurban_pop_change = 0
                        else:
                            distance = np.log((udr_unique[i] - udr_unique[0]))
                            deurban_pop_change = distance / max_distance
                        grid_pop = population_year_1d[loc]
                        total_grid_pop = sum(grid_pop)
                        if i == 0:
                            if sum(population_year_1d)>0:
                                population_year_1d = population_year_1d * (
                                        1 - surplus_urban_pop / sum(population_year_1d))
                            loc = udr_indices >= 0
                            surplus_urban_pop = 0
                        elif (total_grid_pop * deurban_pop_change > surplus_urban_pop):
                            grid_pop_m = grid_pop * (surplus_urban_pop / total_grid_pop)
                            population_year_1d[loc] = grid_pop_m
                            surplus_urban_pop = 0
                        else:
                            new_population = np.rint(grid_pop * (1 - deurban_pop_change))
                            population_year_1d[loc] = new_population
                            surplus_urban_pop = surplus_urban_pop - np.rint(total_grid_pop * deurban_pop_change)

                        # remove urban fabric if new population is below threshold
                        clc_year_1d[loc & (population_year_1d < urban_remove_threshold)] = 51

                elif year >= 2011:
                    # if there should be less population than there is, reduce population by the same proportion
                    decrease_ratio = urban_pop_year / sum(population_year_1d)
                    population_year_1d = population_year_1d * decrease_ratio

            elif surplus_urban_pop < 0:
                if year <= 2011:
                    # if there should be more population than there is, increase population by the same proportion
                    increase_ratio = urban_pop_year / sum(population_year_1d)
                    population_year_1d = population_year_1d * increase_ratio
                else:
                    # add urban areas in available rural space (non-artificial & habitable)
                    rural_mask_lau_empty = (local_pop_density == vLAU) & (((clc_year >= 12) &
                                           (clc_year <= 29)) | (clc_year == 9) | (
                                           clc_year == 32)) & (airports >= year) & (reservoirs >= year)
                    if sum(sum(
                            rural_mask_lau_empty)) < 2:  # case where there is no empty space or only one cell for build-up
                        increase_ratio = urban_pop_year / sum(population_year_1d)
                        population_year_1d = population_year_1d * increase_ratio
                    else:
                        clc_year_1d_r = clc_year[rural_mask_lau_empty].flatten()
                        population_year_1d_r = population_year[rural_mask_lau_empty].flatten()
                        bnp_unique, bnp_indices = np.unique(landuse_bn_to_urban[rural_mask_lau_empty],
                                                            return_inverse=True)
                        cud_r = combined_urban_distance[rural_mask_lau_empty].flatten()
                        i = bnp_unique.size
                        max_density = population_year_1d.max()
                        if cud_r.max() - cud_r.min() <= 1:
                            max_distance = 1
                        else:
                            max_distance = np.log((cud_r.max() - cud_r.min()))
                        while surplus_urban_pop < 0:
                            i = i - 1
                            if vLAU==36:
                                a=1
                            if i < 0:
                                locs = bnp_indices >= i
                                population_year_1d_r[locs] = population_year_1d_r[locs] + (
                                        -surplus_urban_pop / sum(locs))
                                surplus_urban_pop = 0
                            else:
                                loc = bnp_indices == i
                                cud_r_i = cud_r[loc]
                                cud_ranks = rankdata(cud_r_i, method='ordinal')
                                clc_year_1d_r_i = clc_year_1d_r[loc]
                                population_year_1d_r_i = population_year_1d_r[loc]
                                for d in range(1,cud_ranks.size):
                                    loc2 = cud_ranks == d
                                    rank_value = cud_r_i[loc2]
                                    if (rank_value - cud_r.min()) <= 1:
                                        distance = 0
                                    else:
                                        distance = np.log(rank_value - cud_r.min())
                                    urban_pop_addition = (1 - distance / max_distance) * max_density
                                    if -surplus_urban_pop < urban_pop_addition:
                                        population_year_1d_r_i[loc2] = -surplus_urban_pop
                                        surplus_urban_pop = 0
                                        break
                                    else:
                                        existing_pop = population_year_1d_r_i[loc2]
                                        if urban_pop_addition > existing_pop:
                                            population_year_1d_r_i[loc2] = urban_pop_addition
                                            if urban_pop_addition > urban_add_threshold:
                                                clc_year_1d_r_i[loc2] = 2
                                        surplus_urban_pop = surplus_urban_pop + urban_pop_addition - existing_pop

                                population_year_1d_r[loc] = population_year_1d_r_i
                                clc_year_1d_r[loc] = clc_year_1d_r_i

                        population_year[rural_mask_lau_empty] = population_year_1d_r
                        clc_year[rural_mask_lau_empty] = clc_year_1d_r
            # save data
            population_year[urban_mask_lau] = population_year_1d
            clc_year[urban_mask_lau] = clc_year_1d

    return population_year, clc_year

# Rural population
def model_rural_population(clc_year, population_year, population_baseline, vLAUs, vLAU_PG, local_pop_density,
                           pop_baseline_mask):

    clc_mask = (((clc_year >= 3) & (clc_year <= 50)) | ((clc_year >= 52) & (clc_year <= 55))) & pop_baseline_mask
    for vLAU in vLAUs:
        idx = vLAU_PG[:,0] == vLAU
        if vLAU_PG[idx ,[7]] > 0: # check if there is any rural population
            rural_mask_lau_r = (local_pop_density == vLAU) & clc_mask
            # change population by the same proportion
            rural_pop_baseline_r = population_baseline[rural_mask_lau_r].flatten()
            rural_pop_year = np.rint(vLAU_PG[idx, [7]])
            rural_pop_change = rural_pop_year / sum(rural_pop_baseline_r)
            rural_pop_year_r = np.floor(rural_pop_baseline_r * rural_pop_change)
            rural_pop_diff = rural_pop_year - sum(rural_pop_year_r)
            if rural_pop_diff > 0:
                rural_pop_mod = rankdata(np.mod(rural_pop_baseline_r * rural_pop_change, 1), method='ordinal')
                ix = rural_pop_mod <= rural_pop_diff
                rural_pop_year_r[ix] = rural_pop_year_r[ix] + 1
            # save data
            population_year[rural_mask_lau_r] = rural_pop_year_r

    return population_year