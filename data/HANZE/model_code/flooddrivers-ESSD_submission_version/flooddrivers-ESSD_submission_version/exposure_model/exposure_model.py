import copy, re
import geopandas
import numpy as np
import os, sys, inspect
import pandas as pd
import rasterio
from land_use_models import model_airports, model_reservoirs, model_industry, model_infrastructure,  \
    model_construction, model_greenurban, model_burnt_areas, model_agriculture, model_natural_areas_base, \
    model_natural_areas, model_soil_sealing, model_special_cases
from population_models import model_suregional_population, model_urban_population, model_rural_population
from economic_disaggregation import disaggregate_economic

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file
from auxiliary_functions import load_dataset_by_NUTS, prepare_copulas_for_population, NUTS_mask, save_raster_data, \
    write_empty_raster, flood_damage_estimation

def main_model(years, regions_option, hazard_mask_type, uncertainty_option, compute_damage):

    # COMMON PREPROCESSING

    # Load NUTS regions
    NUTS2010 = geopandas.read_file(file('NUTS2010_shapefile'))
    NUTS_code = NUTS2010['Code']

    # initialize random generator and create optional uncertainty sample
    rng = np.random.default_rng(12345)

    # prepare copulas
    LAU_population_data = np.genfromtxt(file('Population_LAU'), delimiter=',')
    copula_samples = rng.uniform(low=0, high=1, size=10000)
    copula_rural = prepare_copulas_for_population(LAU_population_data, 0, 2000, 'normal',5)
    copula_city = prepare_copulas_for_population(LAU_population_data, 1500, 100000, 'Frank',9)

    # define raster datasets
    nuts_dataset = rasterio.open(file('NUTS2010_raster'))
    pop_dataset = rasterio.open(file('Population_100m_baseline'))
    clc_dataset = rasterio.open(file('Land_cover_baseline'))
    imp_dataset = rasterio.open(file('Imperviousness2012'))
    distance_dataset1 = rasterio.open(file('Urban_distance_UN_agglomerations'))
    distance_dataset2 = rasterio.open(file('Urban_distance_UrbanAudit'))
    distance_dataset3 = rasterio.open(file('Urban_distance_Clusters_high_density'))
    distance_dataset5 = rasterio.open(file('Urban_distance_CLC2012'))
    distance_dataset6 = rasterio.open(file('Kernel_density_population'))
    virtual_laus_pop_den = rasterio.open(file('Population_density_vLAU'))
    airport_dataset = rasterio.open(file('Airports'))
    reservoir_dataset = rasterio.open(file('Reservoirs'))
    greenurban_dataset = rasterio.open(file('Selected_green_urban_patches'))
    industry_dataset = rasterio.open(file('Industrial_patches_centroids'))
    landuse_bn_dataset1 = rasterio.open(file('BN_to_urban'))
    landuse_bn_dataset2 = rasterio.open(file('BN_to_crop'))
    landuse_bn_dataset3 = rasterio.open(file('BN_to_past'))
    landuse_bn_dataset4 = rasterio.open(file('BN_from_crop'))
    landuse_bn_dataset5 = rasterio.open(file('BN_from_past'))
    landuse_bn_dataset6 = rasterio.open(file('BN_to_crop_p'))
    landuse_bn_dataset7 = rasterio.open(file('BN_to_past_p'))
    landuse_bn_dataset8 = rasterio.open(file('BN_from_crop_p'))
    landuse_bn_dataset9 = rasterio.open(file('BN_from_past_p'))
    slope_dataset = rasterio.open(file('Slope'))
    polder_dataset = rasterio.open(file('Dutch_polders'))
    street_dataset = rasterio.open(file('ESM2012_streets'))

    # load tabular data
    pop_excel_data = file('NUTS3_database_population_landuse')
    eco_excel_data = file('NUTS3_database_economy')
    NUTS3_population = pd.read_excel(open(pop_excel_data, 'rb'),sheet_name='Population',index_col='Code')
    NUTS3_urban_fraction = pd.read_excel(open(pop_excel_data, 'rb'),sheet_name='Urban fraction',index_col='Code')
    NUTS3_pph = pd.read_excel(open(pop_excel_data, 'rb'), sheet_name='Persons per household', index_col='Code')
    NUTS3_cropland_prc = pd.read_excel(open(pop_excel_data, 'rb'),sheet_name='Croplands',index_col='Code')
    NUTS3_pasture_prc = pd.read_excel(open(pop_excel_data, 'rb'),sheet_name='Pastures',index_col='Code')
    NUTS3_forest_prc = pd.read_excel(open(pop_excel_data, 'rb'),sheet_name='Forests',index_col='Code')
    NUTS3_infrastructure = pd.read_excel(open(pop_excel_data, 'rb'),sheet_name='Infrastructure',index_col='Code')
    NUTS3_GDP = pd.read_excel(open(eco_excel_data, 'rb'), sheet_name='GDP (mln of 2020 euro)', index_col='Code')
    NUTS3_agr_share = pd.read_excel(open(eco_excel_data, 'rb'), sheet_name='% GDP from agriculture', index_col='Code')
    NUTS3_ind_share = pd.read_excel(open(eco_excel_data, 'rb'), sheet_name='% GDP from industry', index_col='Code')
    Fixed_assets_all = pd.read_excel(open(eco_excel_data, 'rb'), sheet_name='Fixed assets', index_col='Code')
    Sector_indices = pd.read_excel(open(eco_excel_data, 'rb'), sheet_name='Sector indices', index_col='Code')

    # NUTS3 region to analyse
    NUTS3_regions = pd.DataFrame(NUTS_code)
    if isinstance(regions_option, str):
        NUTS3_regions_s = NUTS3_regions.loc[NUTS3_regions['Code'].str.contains(regions_option)]
        if len(NUTS3_regions) > 0:
            NUTS3_regions = NUTS3_regions_s
    elif isinstance(regions_option, int):
        NUTS3_regions = NUTS3_regions.loc[NUTS3_regions.index == regions_option]
    else:
        NUTS3_regions = NUTS3_regions.loc[(NUTS3_regions.index>=regions_option.start)&(NUTS3_regions.index<=regions_option.stop)]

    # calibration factors
    urban_remove_threshold = 9
    urban_add_threshold = 81
    industry_elasticity = 0.45

    ## PREPROCESSING THAT VARIES DEPENDING ON THE OPTION

    # optional mask for the extent of natural disaster/hazard map, for which to compute exposure with uncertainty
    hazard_flag = 0
    uncertainty_flag = 0
    iterations = 1

    if isinstance(hazard_mask_type,str):
        hazard_mask_dataset = rasterio.open(file(hazard_mask_type))
        hazard_flag = 1
        # if there is a mask, check if to use uncertainty
        if (uncertainty_option > 0):
            uncertainty_sample = np.rint(rng.uniform(low=0, high=0.9999, size=uncertainty_option) * 10000)
            uncertainty_flag = 1
            iterations = uncertainty_option

    R = np.zeros(iterations)

    if hazard_flag==0:
        # prepare empty rasters (only when not in hazard or uncertainty mode)
        dims = [clc_dataset.height, clc_dataset.width]
        gdp_fa_profile = nuts_dataset.profile
        gdp_fa_profile.data['dtype'] = 'uint32'
        for year in years:
            write_empty_raster(year, pop_dataset.profile, file('Population_100m_year'), np.uint16, dims)
            write_empty_raster(year, clc_dataset.profile, file('Land_cover_year'), np.int8, dims)
            write_empty_raster(year, imp_dataset.profile, file('Imperviousness_year'), np.int8, dims)
            write_empty_raster(year, gdp_fa_profile, file('GDP_year'), np.uint32, dims)
            write_empty_raster(year, gdp_fa_profile, file('Wealth_year'), np.uint32, dims)

    ## CALCULATION LOOP PER REGION
    for region in NUTS3_regions.index:

        ## COMMON PREPROCESSING PER REGION

        # NUTS region codes
        region_code = NUTS_code[region]
        # find grid cells specific for the NUTS region
        region_mask, location = NUTS_mask(NUTS2010, region, nuts_dataset)
        # load other raster data
        population_baseline = load_dataset_by_NUTS(pop_dataset, location, region_mask, 0)
        clc_baseline = load_dataset_by_NUTS(clc_dataset, location, region_mask, 0)
        urban_distance_1 = load_dataset_by_NUTS(distance_dataset1, location, region_mask, 1)
        urban_distance_2 = load_dataset_by_NUTS(distance_dataset2, location, region_mask, 1)
        urban_distance_3 = load_dataset_by_NUTS(distance_dataset3, location, region_mask, 1)
        urban_distance_5 = load_dataset_by_NUTS(distance_dataset5, location, region_mask, 1)
        urban_distance_6 = load_dataset_by_NUTS(distance_dataset6, location, region_mask, 1)
        local_pop_density = load_dataset_by_NUTS(virtual_laus_pop_den, location, region_mask, 0)
        vLAUs = np.unique(local_pop_density)
        # calculate urban distance and ranks
        combined_urban_distance = (urban_distance_1 * 1 + urban_distance_2 * 2 +
                                   urban_distance_3 * 1.5 + urban_distance_5 * 0.5)/5
        # persons per household
        pph_baseline = NUTS3_pph['2011 C'][region_code]
        # urban population baseline and total
        pop_baseline_mask = population_baseline > 0
        urban_mask = (clc_baseline == 1) | (clc_baseline == 2)
        urban_pop_baseline = NUTS3_urban_fraction['2011 C'][region_code] / 100 * NUTS3_population['2011 C'][region_code] * 1000
        urban_pop_baseline_actual = sum(population_baseline[urban_mask])
        urban_pop_correction = urban_pop_baseline_actual / urban_pop_baseline
        # rural population baseline and total
        rural_pop_baseline = (1 - NUTS3_urban_fraction['2011 C'][region_code] / 100) * NUTS3_population['2011 C'][region_code] * 1000
        rural_pop_baseline_actual = sum(population_baseline[~urban_mask])
        rural_pop_correction = rural_pop_baseline_actual / rural_pop_baseline
        # data for modelling other CLC classes
        industrial_GVA_baseline = NUTS3_GDP[2011][region_code] * (100 - NUTS3_agr_share[2011][region_code])/100 # industry + services
        industrial_distance = load_dataset_by_NUTS(industry_dataset, location, region_mask, 0)
        industrial_mask = clc_baseline == 3
        industrial_area_baseline = sum(sum(industrial_mask))
        infrastructure_mask = clc_baseline == 4
        infrastructure_baseline = sum(sum(infrastructure_mask))
        airports = load_dataset_by_NUTS(airport_dataset, location, region_mask, 0)
        reservoirs = load_dataset_by_NUTS(reservoir_dataset, location, region_mask, 0)
        construction_mask = clc_baseline == 9
        greenurban = load_dataset_by_NUTS(greenurban_dataset, location, region_mask, 0)
        # data for agricultural areas
        landuse_bn_other_to_urban = load_dataset_by_NUTS(landuse_bn_dataset1, location, region_mask, 0)
        landuse_bn_other_to_crop = load_dataset_by_NUTS(landuse_bn_dataset2, location, region_mask, 0)
        landuse_bn_other_to_past = load_dataset_by_NUTS(landuse_bn_dataset3, location, region_mask, 0)
        landuse_bn_crop_to_other = load_dataset_by_NUTS(landuse_bn_dataset4, location, region_mask, 0)
        landuse_bn_past_to_other = load_dataset_by_NUTS(landuse_bn_dataset5, location, region_mask, 0)
        landuse_bn_other_from_crop = load_dataset_by_NUTS(landuse_bn_dataset6, location, region_mask, 0)
        landuse_bn_other_from_past = load_dataset_by_NUTS(landuse_bn_dataset7, location, region_mask, 0)
        landuse_bn_crop_from_other = load_dataset_by_NUTS(landuse_bn_dataset8, location, region_mask, 0)
        landuse_bn_past_from_other = load_dataset_by_NUTS(landuse_bn_dataset9, location, region_mask, 0)
        total_region_area = sum(sum(region_mask))
        slope = load_dataset_by_NUTS(slope_dataset, location, region_mask, 0)
        # data for allocating natural areas
        burnt_mask = clc_baseline == 33
        dominant_natural_classes = model_natural_areas_base(clc_baseline, vLAUs, local_pop_density)
        polders = load_dataset_by_NUTS(polder_dataset, location, region_mask, 0)
        mod_flag = 0
        # soil sealing
        soil_sealing_baseline = load_dataset_by_NUTS(imp_dataset, location, region_mask, 0)
        street_baseline = load_dataset_by_NUTS(street_dataset, location, region_mask, 0)

        ## PREPROCESSING THAT VARIES DEPENDING ON THE OPTION

        # list all flood events
        if compute_damage == 1:
            if hazard_mask_type == 'river_flood_mask':
                Events_model = pd.read_csv(file('NUTS3_River_flood_events_list'))
                wd_files = Events_model.loc[Events_model['NUTS3'] == region_code,]
            elif hazard_mask_type == 'coastal_flood_mask':
                Events_model = pd.read_csv(file('NUTS3_Coastal_flood_events_list'))
                wd_files = Events_model.loc[Events_model['NUTS3'] == region_code,]
            elif hazard_mask_type == 'coastal_flood_mask_cf':
                Events_model = pd.read_csv(file('NUTS3_Coastal_flood_events_list_cf'))
                wd_files = Events_model.loc[Events_model['NUTS3'] == region_code,]
            else:
                wd_files = list()

            if len(wd_files)==0:
                print('No events in region: ' + region_code)
                continue
            else:
                R_Fat = np.zeros([len(wd_files), len(years)])
                R_Pop = np.zeros([len(wd_files), len(years)])
                R_FA = np.zeros([len(wd_files), len(years)])
        # Hazard mask (optional, but required for uncertainty analysis)
        elif hazard_flag == 1:
            hazard_mask = load_dataset_by_NUTS(hazard_mask_dataset, location, region_mask, 0)
            # stop running the function if the region is not covered by the hazard mask
            if hazard_mask.max() == 0:
                print('No hazard in region: '+region_code)
                return
            R_Pop = np.zeros([iterations, len(years)])
            R_GDP = np.zeros([iterations, len(years)])
            R_FA = np.zeros([iterations, len(years)])

        ## OPTIONAL UNCERTAINTY LOOP (ONE ITERATION ONLY IS MADE iF UNCERTAINTY OPTION IS NOT USED)
        for idr, r in enumerate(R):

            # samples for population uncertainty analysis (optional)
            if uncertainty_flag == 1:
                copula_sample = int(uncertainty_sample[idr])
            else:
                copula_sample = 1

            # subregional (vLAU) population growth
            vLAU_PG = model_suregional_population(population_baseline, vLAUs, local_pop_density, urban_mask,
                                                  copula_rural, copula_city, copula_samples, urban_distance_6,
                                                  copula_sample, uncertainty_flag)

            ## CALCULATION LOOP PER EACH YEAR
            for idy, year in enumerate(years):

                ## Data preparation

                population_year = population_baseline.astype('float')
                clc_year = copy.deepcopy(clc_baseline)
                soil_sealing_year = copy.deepcopy(soil_sealing_baseline)

                if mod_flag == 1:
                    local_pop_density = load_dataset_by_NUTS(virtual_laus_pop_den, location, region_mask, 0)
                    vLAUs = np.unique(local_pop_density)
                    urban_mask = (clc_baseline == 1) | (clc_baseline == 2)
                    for idx, vLAU in enumerate(vLAUs):
                        urban_mask_lau = (local_pop_density == vLAU) & urban_mask
                        rural_mask_lau = (local_pop_density == vLAU) & (~urban_mask)
                        vLAU_PG[idx, 2] = sum(population_baseline[urban_mask_lau])
                        vLAU_PG[idx, 3] = sum(population_baseline[rural_mask_lau])

                ## Model population and urban fabric

                # # Special cases
                clc_year, population_year, local_pop_density, vLAUs, vLAU_PG, urban_mask, mod_flag = model_special_cases(
                    clc_year, population_year, year, region_code, local_pop_density, vLAUs, vLAU_PG, urban_mask, polders)

                # Urban population and urban_fabric
                population_year, clc_year = model_urban_population(clc_year, population_year, year, region_code,
                                                                   population_baseline, vLAUs, vLAU_PG, local_pop_density,
                                                                   NUTS3_pph, pph_baseline, NUTS3_urban_fraction,
                                                                   NUTS3_population, combined_urban_distance,
                                                                   urban_pop_correction, rural_pop_correction,
                                                                   urban_mask, airports, reservoirs,
                                                                   urban_remove_threshold, urban_add_threshold,
                                                                   landuse_bn_other_to_urban)

                ## Model other CLC classes

                # Airports
                clc_year = model_airports(clc_year, airports, year, population_year)

                # Water
                clc_year = model_reservoirs(clc_year, reservoirs, year, population_year)

                # Rural population
                population_year = model_rural_population(clc_year, population_year, population_baseline, vLAUs,
                                                         vLAU_PG, local_pop_density, pop_baseline_mask)

                # Industrial and commercial
                clc_year = model_industry(clc_year, year, region_code, NUTS3_GDP, NUTS3_agr_share,
                               industrial_GVA_baseline, industrial_area_baseline, industry_elasticity,
                               industrial_distance, industrial_mask, population_year)

                # Roads/railways
                clc_year = model_infrastructure(clc_year, year, region_code, NUTS3_infrastructure,
                                                infrastructure_baseline, infrastructure_mask, combined_urban_distance,
                                                population_year)

                # Construction sites
                clc_year = model_construction(clc_year, year, construction_mask)

                # Green urban areas, sport facilities
                clc_year = model_greenurban(clc_year, greenurban)

                # Croplands
                clc_year = model_agriculture(clc_year, year, region_code, landuse_bn_crop_from_other,
                                             landuse_bn_crop_to_other, landuse_bn_other_from_crop,
                                             landuse_bn_other_to_crop, NUTS3_cropland_prc, total_region_area, slope, 1,
                                             rng, uncertainty_flag)

                # Pastures
                clc_year = model_agriculture(clc_year, year, region_code, landuse_bn_past_from_other,
                                             landuse_bn_past_to_other, landuse_bn_other_from_past,
                                             landuse_bn_other_to_past, NUTS3_pasture_prc, total_region_area, slope, 2,
                                             rng, uncertainty_flag)

                # Burnt areas
                clc_year = model_burnt_areas(clc_year, year, burnt_mask)

                # Natural areas
                clc_year = model_natural_areas(clc_year, year, region_code, vLAUs, local_pop_density,
                                               dominant_natural_classes, NUTS3_forest_prc, total_region_area, polders)

                # Soil sealing (impreviousness)
                soil_sealing_year = model_soil_sealing(soil_sealing_year, clc_year, clc_baseline)

                # disaggregate economic data
                if compute_damage == 1:
                    fatalities, damage_year, pop_exposed = flood_damage_estimation(year, region_code, clc_year,
                                                                  population_year,
                                                                  soil_sealing_year, street_baseline, NUTS3_GDP,
                                                                  NUTS3_agr_share, NUTS3_ind_share, Fixed_assets_all,
                                                                  Sector_indices, wd_files)
                else:
                    gdp_year, wealth_year = disaggregate_economic(year, region_code, clc_year, population_year,
                                                                  soil_sealing_year, street_baseline, NUTS3_GDP,
                                                                  NUTS3_agr_share, NUTS3_ind_share, Fixed_assets_all,
                                                                  Sector_indices)

                # save data or calculate totals per hazard zone
                if hazard_flag==0:
                    # round to integers
                    population_year_r = np.rint(population_year).astype(int)
                    gdp_year_r = np.rint(gdp_year * 1000000).astype(int)
                    wealth_year_r = np.rint(wealth_year * 1000000).astype(int)
                    # save data to files
                    name_suffix = str(year) + '.tif'
                    save_raster_data(file('Population_100m_year') + name_suffix, location, region_mask, population_year_r)
                    save_raster_data(file('Land_cover_year') + name_suffix, location, region_mask, clc_year)
                    save_raster_data(file('Imperviousness_year') + name_suffix, location, region_mask, soil_sealing_year)
                    save_raster_data(file('GDP_year') + name_suffix, location, region_mask, gdp_year_r)
                    save_raster_data(file('Wealth_year') + name_suffix, location, region_mask, wealth_year_r)
                    print(str(region_code) + '_' + str(year))
                elif compute_damage==0:
                    R_Pop[idr, idy] = sum(population_year[hazard_mask == 1])
                    R_GDP[idr, idy] = sum(gdp_year[hazard_mask == 1])*1000
                    R_FA[idr, idy] = sum(wealth_year[hazard_mask == 1])*1000
                    print(str(region_code) + '_' + str(year) + '_' + str(idr))
                elif compute_damage == 1:
                    R_Fat[:, idy] = fatalities
                    R_Pop[:, idy] = pop_exposed
                    R_FA[:, idy] = damage_year * 1000
                    print(str(region_code) + '_' + str(year) + '_' + str(idr))

        # pop uncertainty analysis
        if (hazard_flag == 1) & (uncertainty_flag == 0) & (compute_damage == 0):
            RR = np.concatenate([R_Pop,R_GDP,R_FA])
            np.savetxt(file('Exposure_hazard_zone') + region_code + '.txt', RR, fmt='%d', delimiter=',', )
        if (uncertainty_flag == 1) & (compute_damage == 0):
            R1 = np.percentile(R_Pop,[5,20,50,80,95], axis=0)
            np.savetxt(file('Population_uncertainty') + region_code + '.txt', R1, fmt='%d', delimiter=',',)
            R2 = np.percentile(R_GDP, [5, 20, 50, 80, 95], axis=0)
            np.savetxt(file('GDP_uncertainty') + region_code + '.txt', R2, fmt='%d', delimiter=',', )
            R3 = np.percentile(R_FA, [5, 20, 50, 80, 95], axis=0)
            np.savetxt(file('Wealth_uncertainty') + region_code + '.txt', R3, fmt='%d', delimiter=',', )
        if (compute_damage == 1):
            dates = np.zeros([len(wd_files), 4])
            event_info = np.zeros([len(wd_files), 2])
            for kf, f in enumerate(wd_files.index):
                seg = re.split('_|\.', wd_files.Path[f])
                ep = len(seg)
                event_info[kf, :] = np.array([int(seg[ep - 7]), wd_files.CountryEventID[f]])
                dates[kf, :] = np.array([int(seg[ep - 6]), int(seg[ep - 5]), int(seg[ep - 3]), int(seg[ep - 2])])
            RF_FA = np.concatenate([event_info, R_FA, dates], axis=1)
            RF_Pop = np.concatenate([event_info, R_Pop, dates], axis=1)
            RF_Fat = np.concatenate([event_info, R_Fat, dates], axis=1)
            if hazard_mask_type == 'river_flood_mask':
                np.savetxt(file('River_flood_events_damage') + region_code + '.txt', RF_FA, fmt='%d', delimiter=',', )
                np.savetxt(file('River_flood_events_pop') + region_code + '.txt', RF_Pop, fmt='%d', delimiter=',', )
                np.savetxt(file('River_flood_events_fatalities') + region_code + '.txt', RF_Fat, fmt='%d', delimiter=',', )
            elif hazard_mask_type == 'coastal_flood_mask':
                np.savetxt(file('Coastal_flood_events_damage') + region_code + '.txt', RF_FA, fmt='%d', delimiter=',', )
                np.savetxt(file('Coastal_flood_events_pop') + region_code + '.txt', RF_Pop, fmt='%d', delimiter=',', )
                np.savetxt(file('Coastal_flood_events_fatalities') + region_code + '.txt', RF_Fat, fmt='%d', delimiter=',', )
            elif hazard_mask_type == 'coastal_flood_mask_cf':
                np.savetxt(file('Coastal_flood_events_damage') + '_counterfactual_' + region_code + '.txt', RF_FA, fmt='%d', delimiter=',', )
                np.savetxt(file('Coastal_flood_events_pop') + '_counterfactual_' + region_code + '.txt', RF_Pop, fmt='%d', delimiter=',', )
                np.savetxt(file('Coastal_flood_events_fatalities') + '_counterfactual_' + region_code + '.txt', RF_Fat, fmt='%d', delimiter=',', )
