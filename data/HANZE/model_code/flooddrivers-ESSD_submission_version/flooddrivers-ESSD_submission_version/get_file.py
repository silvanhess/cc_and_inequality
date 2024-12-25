# function to get files used by the model

## DEFINE MAIN PATH TO DOWNLOADED DATA
repo_path = '/p/projects/flooddrivers/'
main_path = repo_path + 'exposure_model/'
# repo_path = 'C:/'
# main_path = repo_path

# paths
path_exposure_tables = main_path + 'HANZE2_products/Exposure_inputs/'
path_landcover = main_path + 'HANZE2_rawdata/Land_cover/'
path_pop_inputs = path_exposure_tables + 'Pop_grid/'
path_baseline_maps = main_path + 'HANZE2_products/Baseline_maps/'
path_esm = path_landcover + 'ESM2012/'
path_urban_data = path_landcover + 'Urban/'
path_agriculture = main_path + 'HANZE2_products/Land_cover/Agriculture/'
path_bn_maps = main_path + 'HANZE2_rawdata/Agriculture/BN_maps/'
path_output_maps = main_path + 'HANZE2_rawdata/Output_maps/'
path_uncertainty = path_output_maps + 'Uncertainty_c/'
path_validation_data = main_path + 'HANZE2_rawdata/Exposure_inputs/Muni_data/'
path_hazard = main_path + 'HANZE2_rawdata/Hazard/'
path_projected_maps = path_output_maps + 'CLC_3arcsec/'
path_EFAS_base = path_landcover + 'EFAS/'
path_HYDE_base = main_path + 'HANZE2_temp/HYDE/'
path_reservoir = path_landcover + 'Reservoirs/'
path_water_demand = main_path + 'HANZE2_temp/water_demand/'
path_coastal_water_levels = repo_path +'Delft3D/SLR_GIA/CWL/'
path_river_discharge = repo_path + 'river_floods_f/'
path_vulnerability = main_path + 'HANZE2_products/Vulnerability/'
path_events = main_path + 'HANZE2_products/Flood_events/'

### TMP
# path_water_demand = main_path + 'HANZE2_temp/EFAS_wateruse/'
# path_coastal_water_levels = main_path + 'HANZE2_products/Surges/Combined_WL/'
# path_river_discharge = main_path + 'HANZE2_temp/EFAS_catchments/'

## define input data
files = {}
# soil sealing (from Copernicus, with adjusted spatial extent)
files['ESM2012_buildings'] = path_esm + "ESM2012_buildings.tif"
files['ESM2012_streets'] = path_esm + "ESM2012_street_ext_adj.tif"
files['Imperviousness2012'] = path_baseline_maps + "IMD2012_extent_adjusted.tif"

# GEOSTAT population 1 km grid (from Eurostat, with adjusted spatial extent)
files['GEOSTAT_1km_population'] = path_pop_inputs + "GEOSTAT_extent_adjusted.tif"

# Corine Land Cover 2012 (from Copernicus, with adjustments)
files['Land_cover_baseline'] = path_baseline_maps + "CLC_base_HANZE2.tif"

# soil sealing - population statistics from preprocessing scripts
files['ESM_GEOSTAT_statistics'] = path_pop_inputs + "ESM_GEOSTAT_statistics.csv"
files['ESM_Street_GEOSTAT_statistics'] = path_pop_inputs + "ESM_Street_GEOSTAT_statistics.csv"
files['IMP_GEOSTAT_statistics'] = path_pop_inputs + "IMP_GEOSTAT_statistics.csv"

# population thresholds generated with preprocessing scripts
files['Population_1km_CLC_pure_cells'] = path_pop_inputs + "Pure_Population_CLC_cells.csv"
files['Population_thresholds'] = path_pop_inputs + "Population_thresholds.csv"

# LAUs and vLAUs (created in GIS)
files['Population_LAU'] = path_urban_data + 'LAU_data.csv'
files['Population_density_vLAU'] = path_urban_data + 'VirtualLAU_PD_int_new.tif'

# NUTS regions (created in GIS)
files['NUTS2010_shapefile'] = path_baseline_maps + 'NUTS2010_final_100m.shp'
files['NUTS2010_raster'] = path_baseline_maps + 'NUTS2010_100m_c.tif'
files['NUTS2021_shapefile'] = main_path + 'HANZE2_products/Admin/NUTS3_regions_v2021_100m.shp'
files['NUTS_2010_to_2021'] = main_path + 'HANZE2_products/Admin/NUTS_change.csv'

# Urban distance datasets (created in GIS)
files['Urban_distance_UN_agglomerations'] = path_urban_data + 'UN_agglomerations_d_int.tif'
files['Urban_distance_UrbanAudit'] = path_urban_data + 'UrbanAudit2018_d_int.tif'
files['Urban_distance_Clusters_high_density'] = path_urban_data + 'Clusters2011_high_density_d_int.tif'
files['Urban_distance_Clusters_urban'] = path_urban_data + 'Clusters2011_urban_d_int.tif'
files['Urban_distance_CLC2012'] = path_urban_data + 'CLC2012_urban_d_int.tif'
files['Kernel_density_population'] = path_urban_data + 'KernelDensityPop1km_int.tif'

# special land-use rasters (created in GIS)
files['Airports'] = main_path + 'HANZE2_products/Land_cover/Airports/Airports_year_v2.tif'
files['Reservoirs'] = path_reservoir + 'Reservoirs_year_v2.tif'
files['Selected_green_urban_patches'] = path_urban_data + 'CLC_141_142_selected_v2.tif'
files['Industrial_patches_centroids'] = path_landcover + 'Industry_infra/Industry_centroids_int.tif'

# Agricultural data from FAO's GAEZ database (projected)
files['GAEZ_wheat'] = path_agriculture + 'ylHr0_whe_LAEA2.tif'
files['GAEZ_grass'] = path_agriculture + 'gras200a_yld_LAEA2.tif'

# Samples for generating the Bayesian Network, created in ArcGIS
files['CLC_changes_sample_data'] = path_agriculture + 'BN/CLC_changes_samples.csv'
files['CLC_nochanges_sample_data'] = path_agriculture + 'BN/CLC_nochanges_samples.csv'

# Bayesian Network probability maps generated with preprocessing scripts
files['BN_sample_data'] = path_agriculture + 'BN/CLC_BN_data.npy'
files['BN_sample_data_validation'] = path_agriculture + 'BN/CLC_BN_data_validation.npy'
files['BN_validation'] = path_agriculture + 'BN/BN_test.txt'
files['BN_map'] = path_bn_maps + 'BN_'
files['BN_to_urban'] = path_bn_maps + 'BN_to_urban.tif'
files['BN_to_crop'] = path_bn_maps + 'BN_to_crop.tif'
files['BN_to_past'] = path_bn_maps + 'BN_to_past.tif'
files['BN_from_crop'] = path_bn_maps + 'BN_from_crop.tif'
files['BN_from_past'] = path_bn_maps + 'BN_from_past.tif'
files['BN_to_crop_p'] = path_bn_maps + 'BN_to_crop_p.tif'
files['BN_to_past_p'] = path_bn_maps + 'BN_to_past_p.tif'
files['BN_from_crop_p'] = path_bn_maps + 'BN_from_crop_p.tif'
files['BN_from_past_p'] = path_bn_maps + 'BN_from_past_p.tif'

# other rasters (adapted from other sources or created in GIS)
files['Slope'] = main_path + 'HANZE2_rawdata/Elevation/Slope_per_mille_int_masked.tif'
files['river_flood_mask'] = path_hazard + 'JRC_flood_maps/JRC_flood_mask_100.tif'
files['coastal_flood_mask'] = path_hazard + 'Delft3D/RAIN_coastalmap_100y.tif'
files['coastal_flood_mask_cf'] = path_hazard + 'Delft3D/RAIN_coastalmap_100y.tif'
files['Dutch_polders'] = path_landcover + 'Special_cases/NL_polders.tif'

# output files generated in the model
files['Population_100m_baseline'] = path_baseline_maps + 'Population_100m.tif'
files['Land_cover_year'] = path_output_maps + 'CLC/CLC_'
files['Population_100m_year'] = path_output_maps + 'Pop/Population_100m_'
files['Imperviousness_year'] = path_output_maps + 'Imp/Imp_'
files['GDP_year'] = path_output_maps + 'GDP/GDP_100m_'
files['Wealth_year'] = path_output_maps + 'FA/FA_100m_'
files['Population_uncertainty'] = path_uncertainty + 'Population_uncertainty_'
files['GDP_uncertainty'] = path_uncertainty + 'GDP_uncertainty_'
files['Wealth_uncertainty'] = path_uncertainty + 'Wealth_uncertainty_'
files['Exposure_hazard_zone'] = path_uncertainty + 'Exposure_hazard_zone_'

# input tabular data (manual data collection)
files['NUTS3_database_population_landuse'] = path_exposure_tables + 'Region_database_population_lu.xlsx'
files['NUTS3_database_economy'] = path_exposure_tables + 'Region_database_economy.xlsx'

# validation data
files['Austria_validation_results_per_LAU'] = path_output_maps + 'Pop/Pop_validation_Austria_'
files['Austria_validation_results_aggregate'] = path_output_maps + 'Pop/Pop_validation_Austria_All.txt'
files['Austria_population_per_LAU'] = path_validation_data + 'LAU2_Austria_pop.shp'
files['Europe_validation_results_per_LAU'] = path_output_maps + 'Pop/Pop_validation_Europe_'
files['Europe_validation_results_aggregate'] = path_output_maps + 'Pop/Pop_validation_Europe_All.txt'
files['Europe_population_per_LAU'] = path_validation_data + 'LAU2_merged_v2021_pop.shp'
files['NUTS3_regions_from_LAUs'] = path_validation_data + 'NUTS3_from_LAU_aggregation.shp'

# visualization
files['Flood_events'] = path_events + 'HANZEv1/Flood_events_v1.0_list_paper.xlsx'
files['Region_exposure_fig'] = path_events + 'Flood_exposure_'
files['Event_exposure_fig'] = path_events + 'Flood_event_'
files['Coastal_zone_expo_path'] = path_output_maps + 'Uncertainty_c/'
files['River_zone_expo_path'] = path_output_maps + 'Uncertainty/'
files['Normalized_losses'] = path_events + 'Exposure_events.csv'
files['Event_pop_unc_fig'] = path_events + 'Population_flood_event_'

# projected maps and output of EFAS land use maps
files['Land_cover_year_p'] = path_projected_maps + 'CLC_'
files['Imperviousness_year_p'] = path_projected_maps + 'Imp_'
files['EFAS_fracforest'] = path_EFAS_base + 'fracforest_European_01min.nc'
files['EFAS_fracirrigated'] = path_EFAS_base + 'fracirrigated_European_01min.nc'
files['EFAS_fracother'] = path_EFAS_base + 'fracother_European_01min.nc'
files['EFAS_fracrice'] = path_EFAS_base + 'fracrice_European_01min.nc'
files['EFAS_fracsealed'] = path_EFAS_base + 'fracsealed_European_01min.nc'
files['EFAS_fracwater'] = path_EFAS_base + 'fracwater_European_01min.nc'
files['EFAS_landmask'] = path_EFAS_base + 'area.nc'
files['EFAS_fracforest_year'] = path_EFAS_base + 'fracforest_European_01min_'
files['EFAS_fracirrigated_year'] = path_EFAS_base + 'fracirrigated_European_01min_'
files['EFAS_fracother_year'] = path_EFAS_base + 'fracother_European_01min_'
files['EFAS_fracrice_year'] = path_EFAS_base + 'fracrice_European_01min_'
files['EFAS_fracsealed_year'] = path_EFAS_base + 'fracsealed_European_01min_'
files['EFAS_fracwater_year'] = path_EFAS_base + 'fracwater_European_01min_'
files['hyde_landlake_dataset'] = path_HYDE_base + 'landlake.asc'
files['hyde_forest_dataset'] = path_HYDE_base + 'forest_wwf_cr.asc'
files['hyde_reservoir_year_dataset'] = path_reservoir + 'HYDE_reservoirs_year.tif'
files['hyde_reservoir_share_dataset'] = path_reservoir + 'HYDE_reservoirs_share.tif'
files['hyde_lu_year_dataset'] = path_HYDE_base + 'zip/'

# data for water demand computation
files['Population_year_p'] = path_projected_maps + 'Population_100m_'
files['water_path'] = path_water_demand

# flood event data
files['TWL'] = path_coastal_water_levels + 'TWL_stations/TWL_'
files['TWL_event'] = path_coastal_water_levels + 'TWL_events/TWL_'
files['TWL_merged_events'] = path_coastal_water_levels + 'Merged_events/Events_'
files['TWL_RP2'] = path_coastal_water_levels + 'TWL_RP2_MaxNoSurge.csv'
files['Coastal_hazard_zones'] = path_coastal_water_levels + 'GIS/Coastal_hazard_zones_cc.shp'
files['TWL_points_raster'] = path_coastal_water_levels + 'GIS/TWL_coastal_points.tif'
files['Coastal_water_depth'] = path_coastal_water_levels + 'GIS/Coastal_water_depth_RP'
files['Coastal_flood_event'] = path_coastal_water_levels + 'Flood_events/'
files['NUTS3_Coastal_flood_events_r'] = path_coastal_water_levels + 'Flood_events_damage/NUTS3_coastal_flood_events_list_'
files['NUTS3_Coastal_flood_events_list'] = path_coastal_water_levels + 'Flood_events_damage/All_NUTS3_coastal_flood_events_list.csv'
files['NUTS3_Coastal_flood_events_list_cf'] = path_coastal_water_levels + 'Flood_events_damage/All_NUTS3_coastal_flood_events_list_counterfactual.csv'
files['Coastal_flood_events_damage'] = path_coastal_water_levels + 'Flood_events_damage/Flood_events_potential_damage_'
files['Coastal_flood_events_pop'] = path_coastal_water_levels + 'Flood_events_damage/Flood_events_population_exposed_'
files['Coastal_flood_events_fatalities'] = path_coastal_water_levels + 'Flood_events_damage/Flood_events_potential_fatalities_'
files['Coastal_flood_event_country'] = path_coastal_water_levels + 'Flood_events_countries/'
files['SeaNames'] = path_coastal_water_levels + 'GIS/SeaNames.shp'
files['SeaGauges'] = path_coastal_water_levels + 'GIS/Coastal_gauges_NUTS.shp'
files['CoastalGaugeObs'] = path_coastal_water_levels + 'Validation/H_'
files['All_Model_Coastal_Flood_Events'] = path_coastal_water_levels + 'All_Model_Coastal_Flood_Events.csv'
files['All_Model_Coastal_Flood_Events_cf'] = path_coastal_water_levels + 'All_Model_Coastal_Flood_Events_cf.csv'
files['All_Model_Coastal_Flood_Events_QObs'] = path_coastal_water_levels + 'All_Model_Coastal_Flood_Events_HObs.csv'
files['Coast_NUTS'] = path_coastal_water_levels + 'JRC_flood_points_NUTS.csv'
files['Q_EFAS'] = path_river_discharge + 'EFAS_outputs/Q_efas_final_'
files['Q_EFAS_aggregate'] = path_river_discharge + 'EFAS_outputs/EFAS_'
files['Q_EFAS_max'] = path_river_discharge + 'EFAS_outputs/EFAS_1995_2004_max.nc'
files['Q_points_EFAS'] = path_river_discharge + 'EFAS_Q_points_hanze.csv'
files['Q_points_NUTS'] = path_river_discharge + 'EFAS_Q_points_hanze_NUTS.csv'
files['Q_selected'] = path_river_discharge + 'Q_annual/EFAS_Q_selected_'
files['Q2'] = path_river_discharge + 'EFAS_outputs/EFAS_Q2.csv'
files['Q_trends'] = path_river_discharge + 'EFAS_outputs/EFAS_Q_trends.csv'
files['Q_event_p'] = path_river_discharge + 'Q_events/'
files['Q_event_c_y'] = path_river_discharge + 'Q_events_c/'
files['River_flood_zones_100_499'] = path_river_discharge + 'GIS/EFAS_thiessen_100_499_cut_LAEA.shp'
files['River_flood_zones_500_plus'] = path_river_discharge + 'GIS/EFAS_thiessen_500plus_cut_LAEA_Q2.shp'
files['Water_depth_EFAS'] = path_river_discharge + 'Aligned_flood_maps/River_EFAS_'
files['Water_depth_RAIN'] = path_river_discharge + 'Aligned_flood_maps/River_RAIN_'
files['Q_points_EFAS_raster'] = path_river_discharge + 'Aligned_flood_maps/EFAS_Q_points.tif'
files['Q_points_RAIN_raster'] = path_river_discharge + 'Aligned_flood_maps/RAIN_Q_points.tif'
files['River_flood_event'] = path_river_discharge + 'Flood_events/'
files['Compound_flood_event'] = path_river_discharge + 'Flood_events/Compound/'
files['NUTS3_River_flood_events_r'] = path_river_discharge + 'Flood_events_damage/NUTS3_river_flood_events_list_'
files['NUTS3_River_flood_events_list'] = path_river_discharge + 'Flood_events_damage/All_NUTS3_river_flood_events_list.csv'
files['River_flood_events_damage'] = path_river_discharge + 'Flood_events_damage/Flood_events_potential_damage_'
files['River_flood_events_pop'] = path_river_discharge + 'Flood_events_damage/Flood_events_population_exposed_'
files['River_flood_events_fatalities'] = path_river_discharge + 'Flood_events_damage/Flood_events_potential_fatalities_'
files['Flood_event_country'] = path_river_discharge + 'Flood_events_countries/'
files['RiverNames'] = path_river_discharge + 'GIS/River_names_NUTS.shp'
files['RiverGauges'] = path_river_discharge + 'Validation/Europe_daily_combined_2022_v2_WGS84_NUTS.shp'
files['RiverGaugeObs'] = path_river_discharge + 'Validation/Q_'
files['All_Model_Flood_Events'] = path_river_discharge + 'All_Model_Flood_Events.csv'
files['All_Model_Flood_Events_QObs'] = path_river_discharge + 'All_Model_Flood_Events_QObs.csv'
files['All_Model_Historical_Events'] = path_river_discharge + 'All_Model_Historical_Events.csv'
files['Model_Events_classified'] = path_river_discharge + 'Events_classified.csv'

# vulnerability analysis
files['VDem'] = path_vulnerability + 'V-Dem-CY-Core-v13.xlsx'
files['All_Model_Flood_Events_SupportDataYE'] = path_vulnerability + 'All_Model_Flood_Events_SupportData_YE.csv'
files['All_Model_Flood_Events_SupportData1950'] = path_vulnerability + 'All_Model_Flood_Events_SupportData_1950.csv'
files['All_Model_Flood_Events_SupportData2020'] = path_vulnerability + 'All_Model_Flood_Events_SupportData_2020.csv'
files['All_Model_Flood_Events_DamageData'] = path_vulnerability + 'All_Model_Flood_Events_DamageData.csv'
files['All_Model_Flood_Events_EventData'] = path_vulnerability + 'All_Model_Flood_Events_EventData.csv'
files['NUTS3_predictor_data'] = path_vulnerability + 'NUTS3_predictors/NUTS3_predictor_data_'
files['MeanPerLabel'] = path_vulnerability + 'MeanData_'
files['Regions_Missing_In_Model'] = path_vulnerability + 'Regions_Missing_In_Model.csv'

# HANZE events
files['Flood_events_v1'] = path_events + 'HANZEv1/Flood_events_v1.0_list.xlsx'
files['Flood_events_v2'] = path_events + 'HANZEv1/Flood_events_v2a_list.xlsx'
files['Flood_events_v2_B_list'] = path_events + 'HANZEv1/Flood_events_v2b_list.xlsx'
files['HistoricalFloods'] = path_river_discharge + 'Validation/'
files['Flood_events_summary_diff'] = path_events + 'Validation/Regions_v1_v2.csv'
files['Flood_events_summary_v2'] = path_events + 'Validation/Regions_v2.csv'
files['HANZE_SHP_2010'] = path_events + 'HANZEv1/HANZE_floods_regions_2010.shp'
files['HANZE_SHP_2021'] = path_events + 'HANZEv1/HANZE_floods_regions_2021.shp'
files['NUTS2010_shapefile_s'] = path_events + 'HANZEv1/Regions_v2010_simplified.shp'
files['NUTS2021_shapefile_s'] = path_events + 'HANZEv1/Regions_v2021_simplified.shp'

### TMP
# files['TWL_event'] = 'C:/HANZE2_temp/TWL_events/TWL_events/TWL_'
# files['Merged_events'] = 'C:/HANZE2_temp/TWL_events/Merged_events/Events_'

# function proper
def file(filename):

    return files[filename]