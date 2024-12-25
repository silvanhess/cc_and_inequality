from scipy.stats import rankdata
import numpy as np

# Airports (CLC 124)
def model_airports(clc_year, airports, year, population_year):
    if year <= 2011:
        airports_mask = (airports > year) & (airports <= 2011)
        clc_year[airports_mask] = 52
    else:
        airports_mask = (airports > 2011) & (airports <= year)
        clc_year[airports_mask] = 6
        population_year[airports_mask] = 0

    return clc_year

# Water - artificial reservoirs (part of CLC 512)
def model_reservoirs(clc_year, reservoirs, year, population_year):
    if year <= 2011:
        reservoir_mask = (reservoirs > year) & (reservoirs <= 2011)
        clc_year[reservoir_mask] = 55
    else:
        reservoir_mask = (reservoirs > 2011) & (reservoirs <= year)
        clc_year[reservoir_mask] = 41
        population_year[reservoir_mask] = 0

    return clc_year

# Industrial and commercial (CLC 121)
def model_industry(clc_year, year, region_code, NUTS3_GDP, NUTS3_agr_share, industrial_GVA_baseline,
                   industrial_area_baseline, industry_elasticity, industrial_distance, industrial_mask,
                   population_year_r):
    # NUTS 3 data on industrial GVA
    industrial_GVA_year = NUTS3_GDP[year][region_code] * (100 - NUTS3_agr_share[year][region_code]) / 100
    industrial_area_change = (industrial_GVA_year / industrial_GVA_baseline) ** industry_elasticity # for calibration
    if industrial_area_change < -1:
        industrial_area_change = -1
    # Industrial CLC class
    industrial_area_year = np.rint(industrial_area_baseline * industrial_area_change)
    surplus_industrial_area = industrial_area_baseline - industrial_area_year

    # removing/adding industrial CLC cells
    if surplus_industrial_area > 0:
        if year <= 2011:
            idr_unique, idr_indices, idr_counts = np.unique(industrial_distance[industrial_mask],
                                                            return_inverse=True, return_counts=True)
            k = idr_unique.size - 1
            locs = 0
            while locs < surplus_industrial_area:
                locs = sum(idr_counts[k:])
                k = k - 1
            loc = idr_indices >= k
            clc_year_1d = clc_year[industrial_mask].flatten()
            clc_year_1d[loc] = 52
            # save industrial data
            clc_year[industrial_mask] = clc_year_1d

    elif surplus_industrial_area < 0:
        if year > 2011:
            rural_mask_empty = (population_year_r == 0) & (((clc_year > 11) & (clc_year < 30)) | (clc_year == 9) | (
                                clc_year == 32))
            if sum(sum(rural_mask_empty)) >0:
                idr_unique, idr_indices, idr_counts = np.unique(industrial_distance[rural_mask_empty],
                                                                return_inverse=True, return_counts=True)
                k = 1
                locs = 0
                while (locs < -surplus_industrial_area) & (k < len(idr_counts) + 2):
                    locs = sum(idr_counts[1:k])
                    k = k + 1
                loc = (idr_indices < (k - 2)) & (idr_indices > 0)
                clc_year_1d = clc_year[rural_mask_empty].flatten()
                clc_year_1d[loc] = 3
                # save industrial data
                clc_year[rural_mask_empty] = clc_year_1d

    return clc_year

# Road/railway infrastructure (CLC 122)
def model_infrastructure(clc_year, year, region_code, NUTS3_infrastructure, infrastructure_baseline,
                         infrastructure_mask, combined_urban_distance, population_year_r):

    # NUTS 3 data on infrastructure area
    infrastructure_year = NUTS3_infrastructure[year][region_code]
    # Infrastructure CLC class
    surplus_infrastructure_area = infrastructure_baseline - infrastructure_year
    # removing/adding infrastructure CLC cells
    if surplus_infrastructure_area > 0:
        ifr_unique, ifr_indices, ifr_counts = np.unique(combined_urban_distance[infrastructure_mask],
                                                        return_inverse=True, return_counts=True)
        k = ifr_unique.size - 1
        locs = 0
        while locs < surplus_infrastructure_area:
            locs = sum(ifr_counts[k:])
            k = k - 1
        loc = ifr_indices >= k
        clc_year_1d = clc_year[infrastructure_mask].flatten()
        clc_year_1d[loc] = 52
        # save infrastructure data
        clc_year[infrastructure_mask] = clc_year_1d

    elif surplus_infrastructure_area < 0:
        rural_mask_empty = clc_year == 9
        ifr_unique, ifr_indices, ifr_counts = np.unique(combined_urban_distance[rural_mask_empty],
                                                        return_inverse=True, return_counts=True)
        if ifr_unique.size > 0:
            if sum(ifr_counts) < -surplus_infrastructure_area:
                k = ifr_counts.size + 1
            else:
                k = 0
                locs = 0
                while locs < -surplus_infrastructure_area:
                    locs = sum(ifr_counts[0:k])
                    k = k + 1
            loc = ifr_indices < (k - 1)
            clc_year_1d = clc_year[rural_mask_empty].flatten()
            clc_year_1d[loc] = 4
            surplus_infrastructure_area = surplus_infrastructure_area + sum(loc)
            # save infrastructure data
            clc_year[rural_mask_empty] = clc_year_1d
        # if there is not enough construction sites, add in any empty space
        if surplus_infrastructure_area < 0:
            rural_mask_empty = (population_year_r == 0) & (((clc_year > 11) & (clc_year < 30)) | (clc_year == 32))
            if sum(sum(rural_mask_empty)) > 0:
                ifr_unique, ifr_indices, ifr_counts = np.unique(combined_urban_distance[rural_mask_empty],
                                                                return_inverse=True, return_counts=True)
                if ifr_indices.size > -surplus_infrastructure_area:
                    k = 0
                    locs = 0
                    while (locs < -surplus_infrastructure_area):
                        locs = sum(ifr_counts[0:k])
                        k = k + 1
                    loc = ifr_indices < (k - 1)
                    clc_year_1d = clc_year[rural_mask_empty].flatten()
                    clc_year_1d[loc] = 4
                    # save infrastructure data
                    clc_year[rural_mask_empty] = clc_year_1d
                else:
                    # if there is not enough space, fill with infrastructure what is available.
                    clc_year[rural_mask_empty] = 4

    return clc_year

# Construction sites (CLC 133)
def model_construction(clc_year, year, construction_mask):

    if year <= 2004:
        clc_year[construction_mask] = 52

    return clc_year

# Green urban areas and sport facilities (CLC 141-142)
def model_greenurban(clc_year, greenurban):

    greenurban_unique = np.unique(greenurban)
    for greenurban_patch in greenurban_unique:
        if greenurban_patch > 0:
            greenurban_mask = greenurban == greenurban_patch
            lowest_clc_class = clc_year[greenurban_mask].min()
            # if there is no more connection with CLC classes 111-122 or 124, remove CLC patch
            if (lowest_clc_class == 5) | (lowest_clc_class > 6):
                greenurban_mask_narrow = (greenurban == greenurban_patch) & (clc_year > 9) & (clc_year < 12)
                clc_year[greenurban_mask_narrow] = 52

    return clc_year

# Burnt areas (CLC 334)
def model_burnt_areas(clc_year, year, burnt_mask):

    if (year <= 2006) | (year >= 2018) :
        clc_year[burnt_mask] = 55

    return clc_year

# Croplands (CLC 211-224, 241-244) and Pasture (CLC 231)
def model_agriculture(clc_year, year, region_code, landuse_bn_agri_from_other, landuse_bn_agri_to_other,
                      landuse_bn_other_from_agri, landuse_bn_other_to_agri, NUTS3_agriculture_prc, total_region_area,
                      slope, variant, rng, uncertainty_flag):

    # croplands, after previous modifications
    if variant ==1:
        agriculture_mask_year = (((clc_year >= 12) & (clc_year <= 17)) | ((clc_year >= 19) & (clc_year <= 22)))
        transition_class = 53
    # pastures, after previous modifications
    else:
        agriculture_mask_year = clc_year == 18
        transition_class = 54
    # NUTS 3 data on cropland/pasture area
    agriculture_year = NUTS3_agriculture_prc[year][region_code] * total_region_area / 100
    # Cropland/pasture CLC class surpluss
    surplus_agriculture_area = sum(sum(agriculture_mask_year)) - agriculture_year
    # removing/adding cropland/pasture CLC cells
    if surplus_agriculture_area > 0:

        # probabilities for past agricultural expansion
        if year <=2011:
            agri_probability = landuse_bn_agri_from_other #8/9
        # probabilities for future agricultural contraction
        else:
            agri_probability = landuse_bn_agri_to_other #4/5

        # remove cropland/pasture according to BN probability of transition
        clc_year_1d_r = clc_year[agriculture_mask_year].flatten()
        bnp_unique, bnp_indices = np.unique(agri_probability[agriculture_mask_year],
                                            return_inverse=True)

        # for uncertainty testing
        if uncertainty_flag==1:
            sample = rng.uniform(low=0, high=1, size=bnp_unique.size)
            bnp_unique_r = bnp_unique * sample
            bnp_unique_r_rank = rankdata(bnp_unique_r, method='ordinal')

        slope_r = slope[agriculture_mask_year].flatten()
        i = bnp_unique.size
        while surplus_agriculture_area > 0:
            i = i - 1
            if i < 0:
                surplus_agriculture_area = 0
            else:
                if uncertainty_flag==1:
                    loc = bnp_indices == bnp_unique_r_rank[i]
                else:
                    loc = bnp_indices == i
                slope_r_i = slope_r[loc]
                if slope_r_i.size < surplus_agriculture_area:
                    clc_year_1d_r[loc] = transition_class
                else:
                    clc_year_1d_r_i = clc_year_1d_r[loc]
                    slope_ranks = rankdata(slope_r_i, method='ordinal')
                    loc2 = slope_ranks > (slope_ranks.size - surplus_agriculture_area + 1)
                    clc_year_1d_r_i[loc2] = transition_class
                    clc_year_1d_r[loc] = clc_year_1d_r_i
                surplus_agriculture_area = surplus_agriculture_area - slope_r_i.size

        # save cropland/pasture data
        clc_year[agriculture_mask_year] = clc_year_1d_r

    elif surplus_agriculture_area < 0:

        # probabilities for past agricultural contraction
        if year <=2011:
            agri_probability = landuse_bn_other_from_agri #6/7
        # probabilities for future agricultural expansion
        else:
            agri_probability = landuse_bn_other_to_agri #2/3

        # area available to croplands
        if variant == 1:
            agriculture_empty = (((clc_year > 22) & (clc_year < 30)) | (clc_year == 18) | ( clc_year == 32) |
                                 (clc_year == 35) | (clc_year == 36) | ((clc_year > 50) & (clc_year < 56)))
            agriculture_class = 12
        # area available to pastures
        else:
            agriculture_empty = (((clc_year > 22) & (clc_year < 30)) | ( clc_year == 32) |
                                 (clc_year == 35) | (clc_year == 36) | ((clc_year > 50) & (clc_year < 56)))
            agriculture_class = 18
        # add cropland according to BN probability of transition
        clc_year_1d_r = clc_year[agriculture_empty].flatten()
        bnp_unique, bnp_indices = np.unique(agri_probability[agriculture_empty],
                                            return_inverse=True)

        # for uncertainty testing
        if uncertainty_flag==1:
            sample = rng.uniform(low=0, high=1, size=bnp_unique.size)
            bnp_unique_r = bnp_unique * sample
            bnp_unique_r_rank = rankdata(bnp_unique_r, method='ordinal')

        slope_r = slope[agriculture_empty].flatten()
        i = bnp_unique.size
        while surplus_agriculture_area < 0:
            i = i - 1
            if i < 0:
                surplus_agriculture_area = 0
            else:
                if uncertainty_flag==1:
                    loc = bnp_indices == bnp_unique_r_rank[i]
                else:
                    loc = bnp_indices == i
                slope_r_i = slope_r[loc]
                if slope_r_i.size < -surplus_agriculture_area:
                    clc_year_1d_r[loc] = agriculture_class
                else:
                    clc_year_1d_r_i = clc_year_1d_r[loc]
                    slope_ranks = rankdata(slope_r_i, method='ordinal')
                    loc2 = slope_ranks <= -surplus_agriculture_area
                    clc_year_1d_r_i[loc2] = agriculture_class
                    clc_year_1d_r[loc] = clc_year_1d_r_i
                surplus_agriculture_area = surplus_agriculture_area + slope_r_i.size

        # save cropland/pasture data
        clc_year[agriculture_empty] = clc_year_1d_r

    return clc_year

def model_natural_areas_base(clc_baseline, vLAUs, local_pop_density):

    # masks of areas of interest
    natural_mask = ((clc_baseline >= 26) & (clc_baseline <= 29)) | ((clc_baseline >= 35) & (clc_baseline <= 38)
                                                                  | (clc_baseline == 32))
    forest_mask = (clc_baseline >= 23) & (clc_baseline <= 25)
    masks = list([natural_mask, forest_mask])

    # create array to save dominant natural classes
    dominant_natural_classes = np.zeros([vLAUs.size, 4])

    # compute dominant forest and non-forest classes
    for k, mask in enumerate(masks):

        # find dominant unique classes
        class_unique, class_counts = np.unique(clc_baseline[mask], return_counts=True)

        # define most frequent natural class in NUTS region as fallback
        if class_counts.size == 0:
            # if the is no natural areas in the region, assign CLC 312 (Coniferous forest) or
            # CLC 324 (Transitional woodland-shrub)
            dominant_class_nuts = (24 if k==0 else 29)
        else:
            dominant_class_nuts = class_unique[class_counts == class_counts.max()]

        # assign most frequent natural class per vLAU
        for idx, vLAU in enumerate(vLAUs):
            # define most frequent natural class in the vLAU
            mask_lau = (local_pop_density == vLAU) & mask
            uncl_unique, uncl_counts = np.unique(clc_baseline[mask_lau], return_counts=True)
            if uncl_counts.size == 0:
                # if the is no natural areas in vLAU, assign most frequent from the region. Count remains 0
                dominant_class_lau = dominant_class_nuts
            else:
                dominant_class_lau = uncl_unique[uncl_counts == uncl_counts.max()]
                dominant_natural_classes[idx, k+2] = uncl_counts.max()

            dominant_natural_classes[idx, k] = \
                (dominant_class_lau if isinstance(dominant_class_lau, int) else dominant_class_lau[0])

    return dominant_natural_classes

def model_natural_areas(clc_year, year, region_code, vLAUs, local_pop_density, dominant_natural_classes,
                        NUTS3_forest_prc, total_region_area, polders):

    # empty spaces where no CLC class was allocated
    unallocated_mask = (clc_year > 50) & (clc_year < 56)

    # special rule for the artificially-created Dutch polders
    if (((region_code=='NL122') | (region_code=='NL230') | (region_code=='NL321')) & (year < 1975)):
        # for polders before the year of construction, change land cover to water bodies (CLC 512)
        polders_mask = (polders > year) & (polders < 1980)
        clc_year[polders_mask] = 41

    # assign most frequent natural class, either forest or non-forest per vLAU
    for idx, vLAU in enumerate(vLAUs):
        n = (0 if dominant_natural_classes[idx,2] > dominant_natural_classes[idx,3] else 1)
        # define most frequent natural class in the vLAU
        unallocated_mask_lau = (local_pop_density == vLAU) & unallocated_mask
        clc_year[unallocated_mask_lau] = dominant_natural_classes[idx,n]

    # forest mask
    forest_mask_year = (clc_year >= 23) & (clc_year <= 25)
    # NUTS 3 data on cropland/pasture area
    forest_year = NUTS3_forest_prc[year][region_code] * total_region_area / 100
    # Cropland/pasture CLC class surpluss
    surplus_forest_area = sum(sum(forest_mask_year)) - forest_year
    # if there is too much/too little forest, modify what class was allocated starting with the most populated LAU
    # (for forest removal) or least (for forest addition)
    k = vLAUs.size
    if surplus_forest_area > 0:
        while (surplus_forest_area > 0) & (k > 0):
            k = k - 1
            if dominant_natural_classes[k,2] <= dominant_natural_classes[k,3]:
                unallocated_mask_lau = (local_pop_density == vLAUs[k]) & unallocated_mask
                clc_year[unallocated_mask_lau] = dominant_natural_classes[k, 0]
                surplus_forest_area = surplus_forest_area - sum(sum(unallocated_mask_lau))
    else:
        k = -1
        while (surplus_forest_area < 0) & (k < vLAUs.size - 1):
            k = k + 1
            if dominant_natural_classes[k,2] > dominant_natural_classes[k,3]:
                unallocated_mask_lau = (local_pop_density == vLAUs[k]) & unallocated_mask
                clc_year[unallocated_mask_lau] = dominant_natural_classes[k, 1]
                surplus_forest_area = surplus_forest_area + sum(sum(unallocated_mask_lau))

    return clc_year

def model_soil_sealing(soil_sealing_year, clc_year, clc_baseline):

    # transitions in CLC classes between baseline and year
    transitions = ~(clc_year == clc_baseline)

    # list of transitions
    to_urban = ((clc_year == 1) | (clc_year == 2)) & transitions # increase to 28%
    to_industry = (clc_year == 3) & transitions  # increase to 45%
    to_infrastructure = (clc_year == 4) & transitions # increase to 29%
    to_airport = (clc_year == 6) & transitions  # increase to 20%
    to_agriculture_from_artificial = (clc_year >= 12) & (clc_year <= 22) & (clc_baseline <= 11) & transitions # reduce to 1%
    to_natural = (clc_year >= 23) & transitions # reduce to 0%

    # modify soil sealing
    soil_sealing_year[to_urban & (soil_sealing_year < 28)] = 28
    soil_sealing_year[to_industry & (soil_sealing_year < 45)] = 45
    soil_sealing_year[to_infrastructure & (soil_sealing_year < 29)] = 29
    soil_sealing_year[to_airport & (soil_sealing_year < 20)] = 20
    soil_sealing_year[to_agriculture_from_artificial & (soil_sealing_year > 1)] = 1
    soil_sealing_year[to_natural & (soil_sealing_year > 0)] = 0

    return soil_sealing_year

def model_special_cases(clc_year, population_year, year, region_code, local_pop_density, vLAUs_old, vLAU_PG,
                        urban_mask, polders):

    if (((region_code=='NL122') | (region_code=='NL230') | (region_code=='NL321')) & (year < 1975)):

        # for polders before the year of construction, assign special 'unallocated land' code 56
        polders_mask = (polders > year) & (polders < 1980)
        clc_year[polders_mask] = 56

        # define urban areas
        urban_mask = (clc_year == 1) | (clc_year == 2)

        # modify virtual LAUs to remove areas that where not yet built
        local_pop_density[polders_mask] = 0
        vLAUs = np.unique(local_pop_density[local_pop_density > 0])
        population_year[polders_mask] = 0
        for idx, vLAU in enumerate(vLAUs_old):
            if np.isin(vLAU,vLAUs):
                urban_mask_lau = (local_pop_density == vLAU) & urban_mask
                rural_mask_lau = (local_pop_density == vLAU) & (~urban_mask)
                vLAU_PG[vLAU_PG[:,0]==vLAU,2] = sum(population_year[urban_mask_lau])
                vLAU_PG[vLAU_PG[:,0]==vLAU,3] = sum(population_year[rural_mask_lau])
            else:
                vLAU_PG[vLAU_PG[:,0]==vLAU,2:4] = 0

        # flag if the "static" data was modified in this timestep"
        mod_flag = 1
    else:
        mod_flag = 0
        vLAUs = vLAUs_old

    return clc_year, population_year, local_pop_density, vLAUs, vLAU_PG, urban_mask, mod_flag