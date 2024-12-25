import pandas as pd
import re
import os, sys, inspect
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

from get_file import file

## this script saves graphs of population (with uncertainty bounds) and a table with normalized losses per event.
## requires exposure uncertainty per hazard zone text files to be computed or downloaded!

# plot settings
plt.ioff()
clrs = sns.color_palette("husl", 5)
plt.rcParams.update({'font.size': 24})

impact_data = file('Flood_events')
Flood_events = pd.read_excel(open(impact_data, 'rb'), sheet_name='List', index_col='No.')
Conversion = pd.read_excel(open(impact_data, 'rb'), sheet_name='Conversion', index_col='Code')
Deflator = pd.read_excel(open(impact_data, 'rb'), sheet_name='Deflator', index_col='Code')
years = np.concatenate([range(1870, 1950, 10), range(1950, 2000, 5), range(2000, 2021)],axis=0)

loss_events = np.zeros([len(Flood_events), 21])

for kf, f in enumerate(Flood_events.itertuples()):

    # check type of event
    type = f.Type
    if type=='Coastal':
        path = file('Coastal_zone_expo_path')
    elif (type=='River'):
        path = file('River_zone_expo_path')
    else:
        continue

    # adjust economic value to 2020 euros
    year = f.Year
    if f.Losses_nom > 0:
        country = f.Code
        currency = f.Currency
        country_currencies = Conversion.loc[Conversion.index==f.Code,]
        factor = country_currencies.loc[country_currencies['Code1']==f.Currency,'Conversion factor']
        euros_nom = f.Losses_nom / factor.values[0]
        deflator = Deflator.loc[Deflator.index==f.Code,year]
        if deflator.values[0] == 0:
            print('No deflator data for event'+str(f.Index))
        euros_2020 = euros_nom / (deflator.values[0] / 100)
    else:
        euros_2020 = np.nan

    regions = re.split(";",f.Regions)

    pop_event = np.zeros([5, len(regions)])
    pop_2020 = np.zeros([5,len(regions)])
    gdp_event = np.zeros([5, len(regions)])
    gdp_2020 = np.zeros([5,len(regions)])
    fas_event = np.zeros([5, len(regions)])
    fas_2020 = np.zeros([5,len(regions)])
    total_pop_expo = np.zeros([3,39])
    for kr,r in enumerate(regions):
        iy = years==2020
        if os.path.isfile(path + 'Wealth_uncertainty_' + r + '.txt'):
            a=1
        else:
            print('No exposure data for event' + str(f.Index))
            continue
        pop_exposure = np.genfromtxt(path + 'Population_uncertainty_' + r + '.txt', delimiter=',')
        gdp_exposure = np.genfromtxt(path + 'GDP_uncertainty_' + r + '.txt', delimiter=',')
        fas_exposure = np.genfromtxt(path + 'Wealth_uncertainty_' + r + '.txt', delimiter=',')

        for u in range(0,5):
            pop_event[u,kr] = np.interp(year,years,pop_exposure[u,])
            pop_2020[u,kr] = pop_exposure[2,iy]
            gdp_event[u,kr] = np.interp(year,years,gdp_exposure[u,])
            gdp_2020[u,kr] = gdp_exposure[2,iy]
            fas_event[u,kr] = np.interp(year,years,fas_exposure[u,])
            fas_2020[u,kr] = fas_exposure[2,iy]

        total_pop_expo += pop_exposure[[0,2,4],]

    fig = plt.figure()
    plt.plot(years, total_pop_expo[1,]/1000, label="Population", c='r')
    plt.fill_between(years, total_pop_expo[0,]/1000, total_pop_expo[2,]/1000, alpha=0.3, facecolor=clrs[0])
    plt.title(f.Country+" "+str(year))
    plt.xticks(np.array([1870,1920,1970,2020]))
    plt.savefig(file('Event_pop_unc_fig')+str(f.Index)+'.png', dpi=120)
    plt.close(fig)

    for u in range(0, 5):
        loss_events[kf,u*4] = f.Fatalities * sum(pop_2020[u,:]) / sum(pop_event[u,:])
        loss_events[kf,u*4+1] = f.Persons * sum(pop_2020[u,:]) / sum(pop_event[u,:])
        loss_events[kf,u*4+2] = euros_2020 * sum(gdp_2020[u,:]) / sum(gdp_event[u,:])
        loss_events[kf,u*4+3] = euros_2020 * sum(fas_2020[u,:]) / sum(fas_event[u,:])

    loss_events[kf, 20] = euros_2020

cols = ['Fat_05','Aff_05','LGDP_05','LFA_05','Fat_20','Aff_20','LGDP_20','LFA_20',
        'Fat_50','Aff_50','LGDP_50','LFA_50','Fat_80','Aff_80','LGDP_80','LFA_80',
        'Fat_95','Aff_95','LGDP_95','LFA_95','Euros2020']

loss_events_df = pd.DataFrame(data=loss_events,index=Flood_events.index, columns=cols)
events_losses = Flood_events.merge(loss_events_df, left_index=True,right_index=True)
events_losses.to_csv(file('Normalized_losses'), sep=',')