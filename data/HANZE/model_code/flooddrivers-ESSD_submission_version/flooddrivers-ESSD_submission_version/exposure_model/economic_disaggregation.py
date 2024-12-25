import copy
import numpy as np

def disaggregate_economic(year, region_code, clc_year, population_year, imp_year, street_baseline, NUTS3_GDP,
                          NUTS3_agr_share,NUTS3_ind_share, Fixed_assets_all, Sector_indices):

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

    # preallocate rasters
    gva_year = np.zeros(clc_year.shape)
    wealth_year = np.zeros(clc_year.shape)

    ## disaggregate GVA
    for k, mask in enumerate(List_of_masks):
        if LU_data[LU_classes[k]][0] > 0:
            # land use component
            if k <= 4:
                # if there are no impervious surfaces, or agriculture/forestry/mining, distribute evenly to land use class
                if (LU_data[LU_classes[k]][3] == 0) | (k <= 1) | (k==3):
                    gva_year[mask] += LU_data[LU_classes[k]][1] * 0.4 / LU_data[LU_classes[k]][0]
                else:
                    gva_year[mask] += imp_year[mask] * \
                                            (LU_data[LU_classes[k]][1] * 0.4 / LU_data[LU_classes[k]][3])
            # population component
            if k <= 1:
                # if there is no agricultural/forest population, distribute evenly to land use class
                if (LU_data[LU_classes[k]][2] == 0):
                    gva_year[mask] += LU_data[LU_classes[k]][1] * 0.4 / LU_data[LU_classes[k]][0]
                else:
                    gva_year[mask] += population_year[mask] * \
                                            (LU_data[LU_classes[k]][1] * 0.6 / LU_data[LU_classes[k]][2])
    # population component of industry, mining and services
    gva_year += population_year * ((industry_services_GVA_year) * 0.6 / Total_population_year)

    ## disaggregate wealth
    for k, mask in enumerate(List_of_masks):
        if (LU_data[LU_classes[k]][0] > 0) & (LU_data[LU_classes[k]][4] > 0):
            if (LU_data[LU_classes[k]][3] == 0):
                wealth_year[mask] += LU_data[LU_classes[k]][4] / LU_data[LU_classes[k]][0]
            else:
                # for infrastructure, combined soil sealing and streets/roads are used if available
                imp_year_mask = (imp_street_year[mask] if (k==6) & (street_flag==1) else imp_year[mask])
                wealth_year[mask] += imp_year_mask * (LU_data[LU_classes[k]][4] / LU_data[LU_classes[k]][3])

    # distribute housing and consumer durables
    wealth_year += population_year * (Residential_assets / Total_population_year)

    return gva_year, wealth_year