# Load Libraries --------------------------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(haven)
library(labelled)
library(countrycode)
library(tmap)
library(sf) # für shape files


# Import data -----------------------------------------------------------------------------------------------------

# survey data

ESS10 <- read_csv("data/ESS/ESS10/ESS10.csv")
ESS8 <- read_csv("data/ESS/ESS8e02_3/ESS8e02_3.csv")



# flood data

HANZE_events <- read_csv("data/HANZE/v2.1/HANZE_events.csv")
# HANZE_references <- read_csv("data/HANZE/datasets/HANZE_references.csv")

# deeznuts

sf_nuts <- read_sf("data/NUTS/NUTS_RG_20M_2010_3035.shp/NUTS_RG_20M_2010_3035.shp")

class(sf_nuts)
sf_nuts$NUTS_ID

nuts <- 
  sf_nuts |> 
  as_tibble() |> 
  rename(region = NUTS_ID) |> 
  select(region,LEVL_CODE, NAME_LATN)

class(nuts)

# Prepare dataset ---------------------------------------------------------------------------------------------

# prepare a dataset with survey participants from wave eight and ten

ESS <- bind_rows(ESS10, ESS8)

ESS_prepared <-
  left_join(x=ESS, y=nuts) |> 
  select(idno, essround, region, LEVL_CODE, NAME_LATN , agea, gndr,lrscale, impenv) |> 
  drop_na(idno, essround, region, agea, gndr,lrscale, impenv) |>
  filter(
    region != "99999",
    lrscale <= 10,
    impenv <= 6
  ) |> 
  mutate(
    gndr = factor(gndr, levels = c(1, 2), labels = c("Male", "Female"))
  )




ESS_prepared |> 
  group_by(LEVL_CODE) |> 
  summarise(count = n()) # change level 3 regions into level 2 regions




# prepare a dataset of floods from 2016 until 2019 in the areas where people were surveyed

HANZE_events_prepared <- 
  HANZE_events |>
    rename(
      end_date = `End date`,
      region = `Regions affected (v2010)`
    ) |> 
    select(ID, end_date, region) |> 
    filter(
      between(as.Date(end_date, format = "%Y-%m-%d"), as.Date("2016-01-01"), as.Date("2019-12-31")) # maybe use year variable instead
    )


# Plot histogram
ggplot(HANZE_events_prepared, aes(x = end_date)) +
  geom_histogram(binwidth = 30, fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of End Dates",
    x = "End Date",
    y = "Count"
  ) +
  theme_minimal()

HANZE_events_prepared_long <- 
  HANZE_events_prepared |>
      separate_rows(region, sep = ";")

# identify the surveyed people that experienced floods
# put them in the treatment group

treatment_group <-     
  left_join(x = HANZE_events_prepared_long, 
            y = ESS_prepared
  ) |> 
  group_by(across(-region)) |> 
  summarize(
    region = paste(unique(region), collapse = ";"),
    .groups = "drop") |> 
  drop_na(impenv)


treatment_group_before <- treatment_group |> filter(essround == 8)
treatment_group_after <- treatment_group |> filter(essround == 10) 


# identify the surveyed people that haven't experienced floods
# put them in the control group


control_group <- ESS_prepared %>%
  filter(!idno %in% treatment_group$idno)

control_group_before <- control_group |> filter(essround == 8)
control_group_after <- control_group |> filter(essround == 10) 


#plot it on a map

sf_HANZE = read_sf("data/HANZE/datasets/HANZE_floods_regions_2021/HANZE_floods_regions_2021.shp")

sf_HANZE_prepared <- 
  sf_HANZE |>
  filter(
    Year >= 2016 &  Year <= 2016
  ) |> 
  mutate(
    "ISO3" = countrycode(Country, origin = 'country.name', destination = 'iso3c')
  ) |> 
  filter(
    ISO3 %in% c(ESS10_prepared_countries$ISO3)
  )



map_HANZE_prepared <- tm_basemap("OpenStreetMap") +
  tm_shape(sf_HANZE_prepared) + 
  tm_dots(col = "blue", size = 0.1)
tmap_mode("view")
tmap_options(check.and.fix = TRUE)
print(map_HANZE_prepared)




# HANZE old version -----------------------------------------------------------------------------------------------

HANZE_v1 <- read_excel("data/HANZE/V1.0/Events_floods.xlsx")

HANZE_v1_prepared <-
  HANZE_v1 |>
  rename(
    end_date = `End date`,
    country_name = `Country name`
  ) |>
  mutate(
    "ISO3" = countrycode(country_name, origin = 'country.name', destination = 'iso3c')
  ) |>
  filter(
    ISO3 %in% c(ESS10_prepared_countries$ISO3),
    between(Year, 2010, 2014)
  )


# NUTS -----------------------------------------------------------------------------------------------------







sf_nuts_prepared <-
  sf_nuts |>
  filter(LEVL_CODE == "2",
         NUTS_ID %in% ESS10_prepared_grouped$region
  )




# FloodArchive ----------------------------------------------------------------------------------------------------


FloodArchive <- read_excel("data/FloodArchive.xlsx")

FloodArchive_prepared <-
  FloodArchive |>
  filter(
    as.Date(Ended, format = "%Y-%m-%d") >= as.Date("2017-01-01")
  ) |>
  mutate(
    ISO3_main = countrycode(Country, origin = 'country.name', destination = 'iso3c'),
    ISO3_other = countrycode(OtherCountry, origin = 'country.name', destination = 'iso3c')
  ) |>
  filter(
    ISO3_main %in% c(ESS10_prepared_countries$ISO3) | ISO3_other %in% c(ESS10_prepared_countries$ISO3)
  )

sf_FloodArchive_prepared <- st_as_sf(FloodArchive_prepared, coords = c("long", "lat"), crs = 4326)
map_FloodArchive <- tm_basemap("OpenStreetMap") +
  tm_shape(sf_FloodArchive_prepared) +
  tm_dots(col = "blue", size = 0.1)
tmap_mode("view")
print(map_FloodArchive)




# emdat -----------------------------------------------------------------------------------------------------------



emdat <- read_excel("data/emdat.xlsx")


emdat_prepared <-
  emdat |> filter(
    ISO %in% c(ESS10_prepared_countries$ISO3),
    `End Year` >= 2016 & `End Year` <= 2020,
    `Disaster Type` == "Flood"
  )


