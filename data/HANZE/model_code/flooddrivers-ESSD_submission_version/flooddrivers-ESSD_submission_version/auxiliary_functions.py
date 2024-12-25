from rasterio.windows import Window
from copulas.bivariate import Frank, Gumbel
from scipy.stats import rankdata, norm, spearmanr
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
import rasterio
import os.path
from netCDF4 import Dataset
from datetime import datetime
import copy, re
import geopandas as gp

from get_file import file

# helper for loading and masking raster for a NUTS region
def load_dataset_by_NUTS(dataset, location, region_mask, distance_adjust):
    read_dataset = dataset.read(1, window=Window(location[0], location[1], location[2], location[3]))
    if distance_adjust == 1:
        read_dataset[read_dataset < 0] = read_dataset.max()
    read_dataset[~region_mask] = 0
    read_dataset[read_dataset < 0] = 0
    return read_dataset

def prepare_copulas_for_population(LAU_pop,threshold_low,threshold_high,copula_type,X_c):

    LAU_data = LAU_pop[(LAU_pop[:, 7] > 0) & (LAU_pop[:, 8] > 9) & (
            LAU_pop[:, 5] >= threshold_low) & (LAU_pop[:, 5] < threshold_high)]
    LAU_data_ranked = rankdata(LAU_data[:, [X_c, 7]], method='min', axis=0) / (LAU_data.shape[0] + 1)
    if copula_type=='Frank':
        copula = Frank()
        copula.fit(LAU_data_ranked)
        param = copula.theta
    elif copula_type=='Gumbel':
        copula = Gumbel()
        r = spearmanr(LAU_data_ranked)
        if r.correlation < 0:
            LAU_data_ranked[:, [1]] = (LAU_data_ranked[:, [1]] - 1)*-1
            copula.fit(LAU_data_ranked)
            param = copula.theta * -1
        else:
            copula.fit(LAU_data_ranked)
            param = copula.theta
    else:
        # return Spearman's correlation for other copulas
        r = spearmanr(LAU_data_ranked)
        param = r.correlation
    X_margins = np.sort(np.stack([LAU_data_ranked[:, 0], LAU_data[:,X_c]], axis=1),axis=0)
    Y_margins = np.sort(np.stack([LAU_data_ranked[:, 1], LAU_data[:,7]], axis=1),axis=0)
    copula = [param,X_margins,Y_margins]

    return copula

def copula_inference_normal(copula,u,v):

    param = copula[0]
    X = copula[1]
    Y = copula[2]

    V = np.interp(v, X[:,1], X[:,0])

    ro = 2 * np.sin( (np.pi / 6 ) * param )
    a = norm.ppf(u, loc=0, scale=1 ) * ( 1 - ro ** 2 ) ** 0.5
    b = ro*norm.ppf(V, loc=0, scale=1)  #v is conditioning variable
    f = norm.cdf(a + b, 0, 1)
    sp = np.interp(f, Y[:,0], Y[:,1])

    return sp

def copula_inference_Frank(copula,u,v):

    theta = copula[0]
    X = copula[1]
    Y = copula[2]

    V = np.interp(v, X[:,1], X[:,0])

    a = -(1 / theta)
    b = 1 - np.exp( - theta)
    c = u ** -1 - 1
    d = np.exp( -theta * V)
    f = a * np.log( 1 - b / (c * d +1 ))
    f[f>1]=1
    sp = np.interp(f, Y[:,0], Y[:,1])

    return sp

def copula_inference_Gumbel(copula,u,v):

    param = copula[0]
    X = copula[1]
    Y = copula[2]

    V = np.interp(v, X[:, 1], X[:, 0])

    flag = 0
    if param<0:
        param = -param
        flag = 1

    C = np.exp(-((-np.log(u))**param+(-np.log(V))**param)**(1/param))
    h = (C * (1./V) * (-np.log(V))**(param-1)) * ((-np.log(u))**param + (-np.log(V))**param)**((1/param)-1)
    h_unique, h_indices = np.unique(h, return_inverse=True)
    uh = np.sort(np.array([h[h_indices], u[h_indices]]))
    f = np.interp(u,uh[0],uh[1])
    if flag == 1:
        f = 1-f

    sp = np.interp(f, Y[:,0], Y[:,1])

    return sp


def sortrows(x):

    x1 = np.zeros([x.shape[0],2])
    x1[:,0] = np.sort(x[:,0])
    for v in x1[:,0]:
        ix = x[:, 0]==v
        ix1 = x1[:, 0] == v
        x1[ix1, 1] = x[ix,1]

    return x1

def NUTS_mask(NUTS2010, region, nuts_dataset):

    EXT_MIN_X = NUTS2010['EXT_MIN_X'][region]
    EXT_MAX_X = NUTS2010['EXT_MAX_X'][region]
    EXT_MIN_Y = NUTS2010['EXT_MIN_Y'][region]
    EXT_MAX_Y = NUTS2010['EXT_MAX_Y'][region]
    gridcode = NUTS2010['gridcode'][region]

    # find the NUTS region in the raster
    start_grid_x = (EXT_MIN_X - 2636000) / 100
    start_grid_y = 40300 - (EXT_MAX_Y - 1386000) / 100
    extent_x = (EXT_MAX_X - EXT_MIN_X) / 100
    extent_y = (EXT_MAX_Y - EXT_MIN_Y) / 100
    location = [start_grid_x, start_grid_y, extent_x, extent_y]
    nuts_region = nuts_dataset.read(1, window=Window(start_grid_x, start_grid_y, extent_x, extent_y))
    # find grid cells specific for the NUTS region
    region_mask = nuts_region == gridcode

    return region_mask, location

def write_empty_raster(variant, profile, filename, data_type, dimensions):

    if profile.data['driver']=='netCDF':
        mode = 'w+'
        ext = '.nc'
    else:
        mode = 'w'
        ext = '.tif'

    full_filename = filename + str(variant) + ext
    if os.path.isfile(full_filename):
        #os.remove(full_filename)
        print(full_filename + " already exists")
    else:
        empty_data = np.zeros(dimensions, dtype = data_type)
        with rasterio.Env():
            with rasterio.open(filename + str(variant) + ext, mode, **profile) as dst:
                dst.write(empty_data, 1)

def save_raster_data(path_and_name, location, region_mask, raster_dataset):

    raster_dataset_year = rasterio.open(path_and_name)
    read_dataset_year = raster_dataset_year.read(1, window=Window(location[0], location[1], location[2], location[3]))
    read_dataset_year[region_mask] = raster_dataset[region_mask]
    profile = raster_dataset_year.profile
    raster_dataset_year.close()
    with rasterio.open(path_and_name, 'r+', **profile) as dst:
        dst.write(read_dataset_year, window=Window(location[0], location[1], location[2], location[3]), indexes=1)
    raster_dataset_year.close()

def save_efas_data(path_and_name, raster_dataset, template_dataset):
    # raster_dataset_year = rasterio.open(path_and_name)
    # profile = raster_dataset_year.profile
    # raster_dataset_year.close()
    # with rasterio.open(path_and_name, 'w+', **profile) as dst:
    #     dst.write(raster_dataset, indexes=1)
    # raster_dataset_year.close()

    if os.path.isfile(path_and_name):
        os.remove(path_and_name)

    # Load template map
    dset = Dataset(template_dataset)
    template_lat = dset.variables['lat'][:]
    template_lon = dset.variables['lon'][:]

    ncfile = Dataset(path_and_name, 'w', format='NETCDF4')
    ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')

    ncfile.createDimension('lon', len(template_lon))
    ncfile.createDimension('lat', len(template_lat))

    ncfile.createVariable('lon', 'f8', ('lon',))
    ncfile.variables['lon'][:] = template_lon
    ncfile.variables['lon'].units = 'degrees_east'
    ncfile.variables['lon'].long_name = 'longitude'

    ncfile.createVariable('lat', 'f8', ('lat',))
    ncfile.variables['lat'][:] = template_lat
    ncfile.variables['lat'].units = 'degrees_north'
    ncfile.variables['lat'].long_name = 'latitude'

    ncfile.createVariable('Band1', np.double, ('lat', 'lon'), zlib=True, chunksizes=(2970, 4530),
                          fill_value=-999999.0)
    ncfile.variables['Band1'][:] = raster_dataset

    return ncfile

def load_hyde_landuse(hyde_path,year):

    # open HYDE
    uopp_dataset = rasterio.open(hyde_path + 'uopp_' + str(year) + 'AD.asc')
    grazing_dataset = rasterio.open(hyde_path + 'grazing' + str(year) + 'AD.asc')
    cropland_dataset = rasterio.open(hyde_path + 'cropland' + str(year) + 'AD.asc')

    HYDE_year_uopp = uopp_dataset.read(1,window=Window(1857, 213, 906, 594)) / 100
    HYDE_year_grazing = grazing_dataset.read(1,window=Window(1857, 213, 906, 594)) / 100
    HYDE_year_cropland = cropland_dataset.read(1,window=Window(1857, 213, 906, 594)) / 100

    return HYDE_year_uopp, HYDE_year_grazing, HYDE_year_cropland

def flood_damage_estimation(year, region_code, clc_year, population_year, imp_year, street_baseline, NUTS3_GDP,
                          NUTS3_agr_share,NUTS3_ind_share, Fixed_assets_all, Sector_indices, wd_files):

    year = year

    # define country of analysis
    country_code = region_code[0:2]

    # load soil sealing
    imp_year[imp_year > 100] = 0
    imp_year = imp_year.astype('float')
    # create combined soil sealing and street/road dataset for infrastructure disaggregation
    imp_street_year = copy.deepcopy(imp_year)
    urban_industrial_mask = clc_year <= 3
    imp_street_year[urban_industrial_mask] = street_baseline[urban_industrial_mask]

    # load population and GDP data
    Total_population_year = sum(sum(population_year))
    GDP_year = NUTS3_GDP[year][region_code]
    Agr_share_year = NUTS3_agr_share[year][region_code] / 100
    Ind_share_year = NUTS3_ind_share[year][region_code] / 100

    # preallocate lists
    List_of_masks = list(np.zeros([7,1]))
    LU_classes = list(['Agriculture', 'Forestry', 'Industry', 'Mining', 'Services_GVA', 'Services_FA', 'Infrastructure'])
    LU_data = {c: np.zeros([5,1]) for c in LU_classes}

    # land cover masks
    List_of_masks[0] = (clc_year >= 12) & (clc_year <= 22) # Agriculture
    List_of_masks[1] = (clc_year >= 23) & (clc_year <= 25) # Forestry
    List_of_masks[2] = clc_year == 3 # Industry
    List_of_masks[3] = clc_year == 7 # Mining
    List_of_masks[4] = ((clc_year > 0) & (clc_year <= 6)) | ((clc_year >= 9) & (clc_year <= 11)) # services (GDP disaggregation)
    List_of_masks[5] = ((clc_year > 0) & (clc_year <= 3)) | ((clc_year >= 9) & (clc_year <= 11)) # services (wealth disaggregation)
    List_of_masks[6] = ((clc_year > 0) & (clc_year <= 6)) # infrastructure

    # Area, population and imperviousness in each class
    street_flag = (1 if sum(imp_street_year[List_of_masks[6]])>0 else 0)
    for k, mask in enumerate(List_of_masks):
        LU_data[LU_classes[k]][0] = sum(sum(mask))
        LU_data[LU_classes[k]][2] = sum(population_year[mask])
        if (k==6) & (street_flag==1):
            # for infrastructure, combined soil sealing and streets/roads are used
            LU_data[LU_classes[k]][3] = sum(imp_street_year[mask])
        else:
            LU_data[LU_classes[k]][3] = sum(imp_year[mask])

    # Gross value added per sector

    # agriculture and forestry GVA and wealth
    agriculture_forest_GVA_year = GDP_year * Agr_share_year
    if (LU_data['Agriculture'][0] + LU_data['Forestry'][0]) > 0:
        forest_intensity = Sector_indices['Forestry'][country_code]
        forest_intensity_area = forest_intensity * LU_data['Forestry'][0] / 100
        forest_share = forest_intensity_area / (LU_data['Agriculture'][0] + forest_intensity_area)
        # save GVA
        LU_data['Agriculture'][1] = agriculture_forest_GVA_year * (1 - forest_share)
        LU_data['Forestry'][1] = agriculture_forest_GVA_year * forest_share
        # save wealth
        LU_data['Agriculture'][4] = LU_data['Agriculture'][1] * Fixed_assets_all[year][country_code][2]
        LU_data['Forestry'][4] = LU_data['Forestry'][1] * Fixed_assets_all[year][country_code][2]
    # services GVA and wealth
    LU_data['Services_GVA'][1] = GDP_year * (1 - Agr_share_year - Ind_share_year)
    LU_data['Services_FA'][4] = LU_data['Services_GVA'][1] * Fixed_assets_all[year][country_code][4]
    industry_services_GVA_year = GDP_year * (1 - Agr_share_year )
    # industry and mining GVA and wealth
    industrial_mining_GVA_year = GDP_year * Ind_share_year
    if (LU_data['Industry'][0] + LU_data['Mining'][0]) > 0:
        mining_intensity = Sector_indices['Mining'][country_code]
        mining_intensity_area = mining_intensity * LU_data['Mining'][0] / 100
        mining_share = mining_intensity_area / (LU_data['Industry'][0] + mining_intensity_area)
        # save GVA
        LU_data['Industry'][1] = industrial_mining_GVA_year * (1 - mining_share)
        LU_data['Mining'][1] = industrial_mining_GVA_year * mining_share
        # save wealth
        LU_data['Industry'][4] = LU_data['Industry'][1] * Fixed_assets_all[year][country_code][3]
        LU_data['Mining'][4] = LU_data['Mining'][1] * Fixed_assets_all[year][country_code][3]

    # in case there is no agricultural and forest area, distribute agriculture with industry
    if (LU_data['Agriculture'][0]+ LU_data['Forestry'][0])==0:
        LU_data['Industry'][1] += agriculture_forest_GVA_year
        LU_data['Industry'][4] += agriculture_forest_GVA_year * Fixed_assets_all[year][country_code][2]
    # in case there is no industrial and mining area, distribute industry with services
    if (LU_data['Industry'][0] + LU_data['Mining'][0])==0:
        LU_data['Services_GVA'][1] += industrial_mining_GVA_year
        LU_data['Services_FA'][4] += industrial_mining_GVA_year * Fixed_assets_all[year][country_code][3]

    ## other wealth data
    LU_data['Infrastructure'][4] = GDP_year * Fixed_assets_all[year][country_code][5]
    Residential_assets = GDP_year * sum(Fixed_assets_all[year][country_code][0:2])

    # define damage functions
    depths = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6])
    List_of_damage_functions = list(np.zeros([8, 1]))
    List_of_damage_functions[0] = interp1d(depths, np.array([0.00, 0.30, 0.55, 0.65, 0.75, 0.85, 0.95, 1.00, 1.00]))
    List_of_damage_functions[1] = interp1d(depths, np.array([0.00, 0.30, 0.55, 0.65, 0.75, 0.85, 0.95, 1.00, 1.00]))
    List_of_damage_functions[2] = interp1d(depths, np.array([0.00, 0.15, 0.27, 0.40, 0.52, 0.70, 0.85, 1.00, 1.00]))
    List_of_damage_functions[3] = interp1d(depths, np.array([0.00, 0.15, 0.27, 0.40, 0.52, 0.70, 0.85, 1.00, 1.00]))
    List_of_damage_functions[5] = interp1d(depths, np.array([0.00, 0.15, 0.30, 0.45, 0.55, 0.75, 0.90, 1.00, 1.00]))
    List_of_damage_functions[6] = interp1d(depths, np.array([0.00, 0.25, 0.42, 0.55, 0.65, 0.80, 0.90, 1.00, 1.00]))
    List_of_damage_functions[7] = interp1d(depths, np.array([0.00, 0.25, 0.40, 0.50, 0.60, 0.75, 0.85, 0.95, 1.00]))

    fatalities = np.zeros([len(wd_files)])
    pop_exposed = np.zeros([len(wd_files)])
    damage_year = np.zeros([len(wd_files)])

    AllEvents = np.unique(wd_files['CountryEventID'].values)
    kf=-1

    for E in AllEvents:
        ie = wd_files['CountryEventID']==E
        Events = wd_files.loc[ie]
        water_depth = np.zeros([population_year.shape[0],population_year.shape[1]])
        for e in Events.index:
            flood_event_dataset = rasterio.open(Events.Path[e])
            water_depth = np.max(np.array([water_depth,flood_event_dataset.read(1)]),axis=0)
            kf = kf+1

        water_depth[water_depth > 6] = 6
        water_depth[water_depth < 0.1] = 0

        pop_exposed[kf] = sum(population_year[water_depth>0])

        # loss of life function from Boyd et al. 2005, "Further specification of the dose-response relationship
        # for flood fatality estimation."
        fatalities[kf] = np.sum(0.34/(1+np.exp(20.37-6.18*water_depth)) * population_year)

        ## disaggregate wealth
        for k, mask in enumerate(List_of_masks):
            if (LU_data[LU_classes[k]][0] > 0) & (LU_data[LU_classes[k]][4] > 0):
                if (LU_data[LU_classes[k]][3] == 0):
                    damage_year[kf] += sum((LU_data[LU_classes[k]][4] / LU_data[LU_classes[k]][0]) \
                                         * List_of_damage_functions[k](water_depth)[mask])
                else:
                    # for infrastructure, combined soil sealing and streets/roads are used if available
                    imp_year_mask = (imp_street_year[mask] if (k==6) & (street_flag==1) else imp_year[mask])
                    damage_year[kf] += sum(imp_year_mask * (LU_data[LU_classes[k]][4] / LU_data[LU_classes[k]][3]) \
                                         * List_of_damage_functions[k](water_depth)[mask])

        # distribute housing and consumer durables
        damage_year[kf] += sum(sum(population_year * (Residential_assets / Total_population_year) * \
                                   List_of_damage_functions[7](water_depth)))

        # # save a compound event
        # if type_flag[kf] == 1:
        #     water_depth[water_depth < 0.1] = flood_event_dataset.profile['nodata']
        #     event_name = f.replace('_events/' + str(int(seg[ep - 7])), '_events/Compound')
        #     event_name = event_name.replace('River', 'Compound')
        #     save_flood_event(event_name, flood_event_dataset.profile, water_depth)

    return fatalities, damage_year, pop_exposed

def get_locations(zone_geo):

    EXT_MIN_X = np.floor(np.round(min(zone_geo['minx'].values))/100)*100
    if EXT_MIN_X < 2636000:
        EXT_MIN_X = 2636000
    EXT_MAX_X = np.ceil(np.round(max(zone_geo['maxx'].values))/100)*100
    if EXT_MAX_X > 6504600:
        EXT_MAX_X = 6504600
    EXT_MIN_Y = np.floor(np.round(min(zone_geo['miny'].values))/100)*100
    if EXT_MIN_Y < 1386000:
        EXT_MIN_Y = 1386000
    EXT_MAX_Y = np.ceil(np.round(max(zone_geo['maxy'].values))/100)*100
    if EXT_MAX_Y > 5416000:
        EXT_MAX_Y = 5416000

    start_grid_x = (EXT_MIN_X - 2636000) / 100
    start_grid_y = 40300 - (EXT_MAX_Y - 1386000) / 100
    extent_x = (EXT_MAX_X - EXT_MIN_X) / 100
    extent_y = (EXT_MAX_Y - EXT_MIN_Y) / 100
    location = [start_grid_x, start_grid_y, extent_x, extent_y]

    return location

def save_flood_event(path_and_name, profile, water_depth):

    if os.path.isfile(path_and_name):
        os.remove(path_and_name)

    with rasterio.Env():
        with rasterio.open(path_and_name, 'w', **profile) as dst:
            dst.write(water_depth, 1)

def nuts_2010_to_2021(Regions2021, FloodExtent, E_Regs2010, NUTS2010_to_2021):

    NUTS2010 = np.array(re.split(';',E_Regs2010))
    NUTS2010 = NUTS2010[NUTS2010.argsort()]
    ix = np.isin(NUTS2010_to_2021[:,0],NUTS2010)
    NUTS2021 = NUTS2010_to_2021[ix,1]

    iz = NUTS2021=='change'

    if np.sum(iz)>0:

        FloodExtent_NUTS2021 = gp.overlay(Regions2021, FloodExtent, how='intersection', keep_geom_type=False)
        FE_area = FloodExtent_NUTS2021.geometry.area
        Regions2021_E = np.unique(FloodExtent_NUTS2021['NUTS'].values)
        Sarea = np.zeros([len(Regions2021_E)])
        for k, r in enumerate(Regions2021_E):
            Sarea[k] = sum(FE_area[FloodExtent_NUTS2021['NUTS'] == r])

        ia = Sarea > 500000
        NUTS2021 = Regions2021_E[ia]

    NUTS2021_sort = np.sort(NUTS2021)

    return NUTS2021_sort

def GDPpercapita(years, Countries, NUTS3_population, NUTS3_GDP):

    NUTS3 = copy.deepcopy(NUTS3_GDP.index.values)
    for k, i in enumerate(NUTS3):
        if isinstance(NUTS3[k], str):
            NUTS3[k] = NUTS3[k][0:2]

    GDPpc = np.zeros([len(Countries.index), len(years)])
    for k, c in enumerate(Countries.index):
        ix = np.isin(NUTS3, c)
        sum_pop = np.sum(NUTS3_population[years][ix], axis=0)
        sum_GDP = np.sum(NUTS3_GDP[years][ix], axis=0)
        GDPpc[k, :] = sum_GDP / sum_pop

    return GDPpc

def NUTS3_data_extract_preprocess(years, NUTS3_population, NUTS3_GDP, NUTS3_Agr, NUTS3_Ind, Country_wealth,
                                  NUTS3_urban_pop, NUTS3_area, LU_stats, VDem_indicators, NUTS_to_countries, NUTS_code):

    Vars = 25
    column_list_lu = ['NUTS3','Urban','OtherArt','Agri','Forest','Other','Imp']

    for r in NUTS_code:

        SupportData = np.zeros([len(years), Vars])
        print(r)

        for k, y in enumerate(years):

            GDP_Sec = np.zeros([3])
            Wealth = np.zeros([6])
            LandUse = np.zeros([6])
            LU_stats_y = pd.DataFrame(data=LU_stats[k,:,:], columns=column_list_lu)
            LU_stats_y = LU_stats_y.set_index('NUTS3')

            Area_r = NUTS3_area[r]
            Pop_r = NUTS3_population[y][r]
            GDP_r = NUTS3_GDP[y][r]
            c = r[0:2]
            GDP_Sec[0] = GDP_r * NUTS3_Agr[y][r] / 100
            GDP_Sec[1] = GDP_r * NUTS3_Ind[y][r] / 100
            GDP_Sec[2] = GDP_r * (100 - NUTS3_Agr[y][r] - NUTS3_Ind[y][r]) / 100
            Wealth[0] = GDP_r * Country_wealth[y][c][0]
            Wealth[1] = GDP_r * Country_wealth[y][c][1]
            Wealth[2] = GDP_r * NUTS3_Agr[y][r] / 100 * Country_wealth[y][c][2]
            Wealth[3] = GDP_r * NUTS3_Ind[y][r] / 100 * Country_wealth[y][c][3]
            Wealth[4] = GDP_r * (100 - NUTS3_Agr[y][r] - NUTS3_Ind[y][r]) / 100 * Country_wealth[y][c][4]
            Wealth[5] = GDP_r * Country_wealth[y][c][5]
            LandUse[0] = LU_stats_y['Urban'][r]
            LandUse[1] = LU_stats_y['OtherArt'][r]
            LandUse[2] = LU_stats_y['Agri'][r]
            LandUse[3] = LU_stats_y['Forest'][r]
            LandUse[4] = LU_stats_y['Other'][r]
            LandUse[5] = LU_stats_y['Imp'][r] * Area_r
            Urban = Pop_r * NUTS3_urban_pop[y][r] / 100

            SupportData[k, 0] = y
            SupportData[k, 1] = GDP_r
            SupportData[k, 2:5] = GDP_Sec
            SupportData[k, 5:11] = Wealth
            SupportData[k, 11:17] = LandUse
            SupportData[k, 17] = Area_r
            SupportData[k, 18] = Urban
            SupportData[k, 19] = Pop_r

        SupportDataYE = np.zeros([121,Vars])
        for i in range(0,Vars-5):
            InterpData = interp1d(years, SupportData[:, i])
            for k,y in enumerate(range(1900,2021)):
                SupportDataYE[k,i] = InterpData(y)
        # V-Dem data
        for k, y in enumerate(range(1900, 2021)):
            VDem = np.zeros([5])
            VDem_country = NUTS_to_countries[y][r]
            VDem_indicators_y = VDem_indicators.loc[VDem_indicators['year'] == y,]
            VDem[0] = VDem_indicators_y['v2x_polyarchy'][VDem_country]
            VDem[1] = VDem_indicators_y['v2x_libdem'][VDem_country]
            VDem[2] = VDem_indicators_y['v2x_partipdem'][VDem_country]
            VDem[3] = VDem_indicators_y['v2x_delibdem'][VDem_country]
            VDem[4] = VDem_indicators_y['v2x_egaldem'][VDem_country]
            SupportDataYE[k, 20:] = VDem

        column_list_sd = ['Year', 'GDP', 'GDP_Agr', 'GDP_Ind', 'GDP_Ser', 'FA_Res', 'FA_Con', 'FA_Agr', 'FA_Ind',
                          'FA_Ser', 'FA_Inf', 'LU_Urban', 'LU_OtherArt', 'LU_Agri', 'LU_Forest', 'LU_Other',
                          'SoilSeal', 'Area', 'PopUrban', 'Pop',
                          'ElectDem', 'LibDem', 'ParticipDem', 'DelibDem', 'EgalDem']

        SupportDataE_YE = pd.DataFrame(data=SupportDataYE, columns=column_list_sd)

        SupportDataE_YE.to_csv(file('NUTS3_predictor_data') + r + '.csv', sep=',')

def NUTS3_data_extract(Regions_E, y, AllEventID):

    Vars = 24
    SupportData = np.zeros([Vars])

    GDP = 0
    Pop = 0
    GDP_Sec = np.zeros([3])
    Wealth = np.zeros([6])
    LandUse = np.zeros([6])
    VDem = np.zeros([5])
    Urban = 0

    for r in Regions_E:

        Predictors = pd.read_csv(file('NUTS3_predictor_data') + r + '.csv', index_col='Year')

        Pop_r = Predictors['Pop'][y]
        GDP_r = Predictors['GDP'][y]
        GDP += GDP_r
        Pop += Pop_r
        GDP_Sec[0] += Predictors['GDP_Agr'][y]
        GDP_Sec[1] += Predictors['GDP_Ind'][y]
        GDP_Sec[2] += Predictors['GDP_Ser'][y]
        Wealth[0] += Predictors['FA_Res'][y]
        Wealth[1] += Predictors['FA_Con'][y]
        Wealth[2] += Predictors['FA_Agr'][y]
        Wealth[3] += Predictors['FA_Ind'][y]
        Wealth[4] += Predictors['FA_Ser'][y]
        Wealth[5] += Predictors['FA_Inf'][y]
        LandUse[0] += Predictors['LU_Urban'][y]
        LandUse[1] += Predictors['LU_OtherArt'][y]
        LandUse[2] += Predictors['LU_Agri'][y]
        LandUse[3] += Predictors['LU_Forest'][y]
        LandUse[4] += Predictors['LU_Other'][y]
        LandUse[5] += Predictors['SoilSeal'][y]
        Urban += Predictors['PopUrban'][y]
        VDem[0] += Pop_r * Predictors['ElectDem'][y]
        VDem[1] += Pop_r * Predictors['LibDem'][y]
        VDem[2] += Pop_r * Predictors['ParticipDem'][y]
        VDem[3] += Pop_r * Predictors['DelibDem'][y]
        VDem[4] += Pop_r * Predictors['EgalDem'][y]

    SupportData[0] = AllEventID
    SupportData[1] = GDP / Pop * 1000
    SupportData[2:5] = GDP_Sec / GDP * 100
    Total_wealth = np.sum(Wealth)
    for w in range(0, 6):
        SupportData[5+w] = Wealth[w] / Total_wealth * 100
    Total_area = np.sum(LandUse)
    for w in range(0, 6):
        SupportData[11+w] = LandUse[w] / Total_area * 100
    SupportData[17] = Urban / Pop * 100
    SupportData[18] = Pop / Total_area * 10000
    for w in range(0, 5):
        SupportData[19+w] = VDem[w] / Pop

    return SupportData

def Aggregate_HANZE_numbers(Reg_count, Regions_E, E_Year, AllEventID, NUTS3):

    HANZE_events_numbers = np.zeros([9])

    ix = np.isin(NUTS3,Regions_E)
    n_regions = len(Regions_E)

    HANZE_events_numbers[0] = AllEventID
    HANZE_events_numbers[1] = np.sum(Reg_count[ix, E_Year - 1880:E_Year - 1870])
    HANZE_events_numbers[2] = np.sum(Reg_count[ix, E_Year - 1890:E_Year - 1870])
    HANZE_events_numbers[3] = np.sum(Reg_count[ix, E_Year - 1900:E_Year - 1870])
    HANZE_events_numbers[4] = np.sum(Reg_count[ix, E_Year - 1920:E_Year - 1870])
    HANZE_events_numbers[5] = np.sum(Reg_count[ix, E_Year - 1885:E_Year - 1875])
    HANZE_events_numbers[6] = np.sum(Reg_count[ix, E_Year - 1895:E_Year - 1875])
    HANZE_events_numbers[7] = np.sum(Reg_count[ix, E_Year - 1905:E_Year - 1875])
    HANZE_events_numbers[8] = np.sum(Reg_count[ix, E_Year - 1925:E_Year - 1875])
    HANZE_events_numbers[1:] = HANZE_events_numbers[1:] / n_regions

    return HANZE_events_numbers

def Compute_flood_risk(Regions_E, Type_E, E_Year, E, Coastal_ID, years_c, years_r, Events_NUTS3, C, Events_hanze_v2):

    loss = np.zeros([3])
    loss_E = np.zeros([3])
    risk = np.zeros([3])
    risk_E = np.zeros([3])
    Area_RP = np.zeros([4])
    GDP_1950 = 0.001
    Pop_1950 = 0.001
    GDP_E = 0.001
    Pop_E = 0.001

    for k, r in enumerate(Regions_E):
        d_file_c = file('Coastal_flood_events_damage') + r + '.txt'
        p_file_c = file('Coastal_flood_events_pop') + r + '.txt'
        f_file_c = file('Coastal_flood_events_fatalities') + r + '.txt'
        d_file_r = file('River_flood_events_damage') + r + '.txt'
        p_file_r = file('River_flood_events_pop') + r + '.txt'
        f_file_r = file('River_flood_events_fatalities') + r + '.txt'
        risk_d_c = pd.read_csv(d_file_c, delimiter=',', header=None).values if os.path.isfile(d_file_c) else np.zeros([1,47])
        risk_p_c = pd.read_csv(p_file_c, delimiter=',', header=None).values if os.path.isfile(p_file_c) else np.zeros([1,47])
        risk_f_c = pd.read_csv(f_file_c, delimiter=',', header=None).values if os.path.isfile(f_file_c) else np.zeros([1,47])
        risk_d_r = pd.read_csv(d_file_r, delimiter=',', header=None).values if os.path.isfile(d_file_r) else np.zeros([1,7])
        risk_p_r = pd.read_csv(p_file_r, delimiter=',', header=None).values if os.path.isfile(p_file_r) else np.zeros([1,7])
        risk_f_r = pd.read_csv(f_file_r, delimiter=',', header=None).values if os.path.isfile(f_file_r) else np.zeros([1,7])

        iy_comp = risk_d_r[:, 1] == E
        iy_r = risk_d_r[:, 1] == E - 200000
        if (E >= 100000) & (E < 200000):
            iy_c = risk_d_c[:, 1] == E - 100000
        else:
            iy_c = risk_d_c[:, 1] == Coastal_ID - 100000

        if (sum(iy_comp) + sum(iy_r) + sum(iy_c)) > 0:

            Predictors = pd.read_csv(file('NUTS3_predictor_data') + r + '.csv', index_col='Year')
            Pop_E += Predictors['Pop'][E_Year] * 1000
            GDP_E += Predictors['GDP'][E_Year] * 1000
            Pop_1950 += Predictors['Pop'][1950] * 1000
            GDP_1950 += Predictors['GDP'][1950] * 1000

            # interp1d(years_r, SupportData[:, i])(E_Year)
            if Type_E[k]=='Compound':
                ix = risk_d_r[:, 1] < 200000
                loss[0] += np.sum(risk_d_r[ix,:], axis=0)[2] + np.sum(risk_d_c, axis=0)[2]
                loss[1] += np.sum(risk_p_r[ix,:], axis=0)[2] + np.sum(risk_p_c, axis=0)[2]
                loss[2] += np.sum(risk_f_r[ix,:], axis=0)[2] + np.sum(risk_f_c, axis=0)[2]

                ID_E = E
                if (sum(iy_comp)==0) & (sum(iy_r)==0):
                    loss_E[0] += interp1d(years_c, np.max(risk_d_c[iy_c, 2:len(years_c) + 2], axis=0))(E_Year)
                    loss_E[1] += interp1d(years_c, np.max(risk_p_c[iy_c, 2:len(years_c) + 2], axis=0))(E_Year)
                    loss_E[2] += interp1d(years_c, np.max(risk_f_c[iy_c, 2:len(years_c) + 2], axis=0))(E_Year)
                    ID_E = Coastal_ID
                else:
                    if sum(iy_comp)==0:
                        iy = iy_r
                        ID_E = E - 200000
                    else:
                        iy = iy_comp
                    loss_E[0] += interp1d(years_r, np.max(risk_d_r[iy, 2:len(years_r) + 2], axis=0))(E_Year)
                    loss_E[1] += interp1d(years_r, np.max(risk_p_r[iy, 2:len(years_r) + 2], axis=0))(E_Year)
                    loss_E[2] += interp1d(years_r, np.max(risk_f_r[iy, 2:len(years_r) + 2], axis=0))(E_Year)

            elif Type_E[k]=='Coastal':
                if (E >= 100000) & (E < 200000):
                    Coastal_ID = E
                loss[0] += np.sum(risk_d_c, axis=0)[2]
                loss[1] += np.sum(risk_p_c, axis=0)[2]
                loss[2] += np.sum(risk_f_c, axis=0)[2]
                iy = risk_d_c[:, 1] == Coastal_ID - 100000
                loss_E[0] += interp1d(years_c, np.max(risk_d_c[iy, 2:len(years_c) + 2], axis=0))(E_Year)
                loss_E[1] += interp1d(years_c, np.max(risk_p_c[iy, 2:len(years_c) + 2], axis=0))(E_Year)
                loss_E[2] += interp1d(years_c, np.max(risk_f_c[iy, 2:len(years_c) + 2], axis=0))(E_Year)
                ID_E = Coastal_ID
            else:
                ix = risk_d_r[:, 1] < 200000
                loss[0] += np.sum(risk_d_r[ix, :], axis=0)[2]
                loss[1] += np.sum(risk_p_r[ix, :], axis=0)[2]
                loss[2] += np.sum(risk_f_r[ix, :], axis=0)[2]
                if E >= 200000:
                    E = E - 200000
                iy = risk_d_r[:, 1] == E
                loss_E[0] += interp1d(years_r, np.max(risk_d_r[iy, 2:len(years_r) + 2], axis=0))(E_Year)
                loss_E[1] += interp1d(years_r, np.max(risk_p_r[iy, 2:len(years_r) + 2], axis=0))(E_Year)
                loss_E[2] += interp1d(years_r, np.max(risk_f_r[iy, 2:len(years_r) + 2], axis=0))(E_Year)
                ID_E = E

            Events_NUTS3_E0 = Events_NUTS3.loc[(Events_NUTS3['CountryEventID']==ID_E) &
                                              (Events_NUTS3['NUTS3']==r),:]
            E_Area = max(Events_NUTS3_E0['Area'].values)
            Events_NUTS3_E = Events_NUTS3_E0.loc[Events_NUTS3_E0['Area'] == E_Area,]
            Area_RP[0] += E_Area
            Area_RP[1] += max(Events_NUTS3_E['RPv1'].values * Events_NUTS3_E['Area'].values)
            Area_RP[2] += max(Events_NUTS3_E['RPv2'].values * Events_NUTS3_E['Area'].values)
            Area_RP[3] += max(Events_NUTS3_E['WaterDepth'].values * Events_NUTS3_E['Area'].values)

    risk[0] = loss[0] / GDP_1950
    risk[1] = loss[1] / Pop_1950
    risk[2] = loss[2] / Pop_1950

    risk_E[0] = loss_E[0] / GDP_E * 100
    risk_E[1] = loss_E[1] / Pop_E * 100
    risk_E[2] = loss_E[2] / Pop_E * 100

    risk = risk * 100 / 71

    Area_RP[1] = Area_RP[1] / Area_RP[0]
    Area_RP[2] = Area_RP[2] / Area_RP[0]
    Area_RP[3] = Area_RP[3] / Area_RP[0] / 10 # from mm to cm

    loss_E[loss_E==0] = 1 # add to avoid NaNs when computing e.g. relative fatalities, when reported fatalities = 0

    ReportedLosses = np.full([8], np.nan)
    if ((C>0) & (C < 3000)) | ((C>=6000) & (C < 7000)):
        ReportedLosses[0] = Events_hanze_v2['Area'][int(C)] * 100
        ReportedLosses[1] = Events_hanze_v2['Fatalities'][int(C)]
        ReportedLosses[2] = Events_hanze_v2['Persons'][int(C)]
        ReportedLosses[3] = Events_hanze_v2['Losses_norm'][int(C)] * 1000
        ReportedLosses[4] = ReportedLosses[0] / Area_RP[0] * 100
        ReportedLosses[5] = ReportedLosses[1] / loss_E[2] * 100
        ReportedLosses[6] = ReportedLosses[2] / loss_E[1] * 100
        ReportedLosses[7] = ReportedLosses[3] / loss_E[0] * 100

    risk_data = np.concatenate([risk, risk_E, loss_E, ReportedLosses, Area_RP])

    return risk_data

def con_strings(string, option):
    E_String = ''
    for r in string:
        if isinstance(option,int):
            E_String = E_String + str(r) + ';'
        else:
            E_String = E_String + option.loc[option['gridcode'] == r, 'Code'].values[0] + ';'
    E_String2 = E_String[0:len(E_String) - 1]

    return E_String2