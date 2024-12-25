# Load Libraries --------------------------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(haven)
library(labelled)
library(countrycode)
library(tmap)
library(sf) # für shape files

# NUTS -----------------------------------------------------------------------------------------------------




nuts <- read_excel("data/nuts.xlsx", sheet = "NUTS2024")
sf_nuts <- read_sf("data/NUTS_RG_20M_2010_3035.shp/NUTS_RG_20M_2010_3035.shp")


sf_nuts_prepared <- 
  sf_nuts |> 
  filter(LEVL_CODE == "2",
         NUTS_ID %in% ESS10_prepared_grouped$region
  )



# Prepare Survey Data ---------------------------------------------------------------------------------------------

ESS10 <- read_csv("data/ESS10.csv")

ESS10_prepared <- 
  ESS10 |>
  left_join(x=ESS10, y=nuts, by = c(region = "NUTS Code")) |> 
  select(cntry, region, agea, gndr, impenv) |> 
  drop_na() |> 
  filter(
    region != "99999",
    impenv <= 6
  ) |> 
  mutate(
    "ISO3" = countrycode(cntry, origin = "iso2c", destination = "iso3c"),
    gndr = factor(gndr, levels = c(1, 2), labels = c("Male", "Female"))
  )

ESS10_prepared_grouped <- 
  ESS10_prepared |> 
    group_by(region) |> 
    summarise(
      impenv_mean = mean(impenv)
    )

ESS10_prepared_countries <- 
  ESS10_prepared |> 
  group_by(ISO3) |> 
  summarise(
    count = n()
  )


# countries that were surveyed in both wave 10 and 11
# ESS_countries <- inner_join(
#   x = ESS11_countries, 
#   y = ESS10_countries, 
#   by = "cntry") |> 
#   mutate(
#     ISO3 = countrycode(cntry, origin = "iso2c", destination = "iso3c")
#   )


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




# HANZE-------------------------------------------------------------------------------------------

HANZE_events <- read_csv("data/HANZE/datasets/HANZE_events.csv")
HANZE_references <- read_csv("data/HANZE/datasets/HANZE_references.csv")
sf_HANZE = read_sf("data/HANZE/datasets/HANZE_floods_regions_2021/HANZE_floods_regions_2021.shp")


HANZE_events_prepared <- 
  HANZE_events |>
    select(-`Regions affected (v2021)`) |> 
    rename(
      end_date = `End date`,
      country_name = `Country name`,
      region = `Regions affected (v2010)`
    ) |> 
    mutate(
      "ISO3" = countrycode(country_name, origin = 'country.name', destination = 'iso3c')
    ) |> 
    filter(
      ISO3 %in% c(ESS10_prepared_countries$ISO3),
      between(as.Date(end_date, format = "%Y-%m-%d"), as.Date("2016-01-01"), as.Date("2019-12-31"))
    )

HANZE_events_prepared_long <- 
  HANZE_events_prepared |>
      separate_rows(region, sep = ";")


HANZE_and_ESS10 <-     
  left_join(x = HANZE_events_prepared_long, 
            y = ESS10_prepared_grouped
  ) |> 
  group_by(across(-region)) |> 
  summarize(
    region = paste(unique(region), collapse = ";"),
    .groups = "drop")


treatment_group <- 
  HANZE_and_ESS10 |> 
    filter(
      impenv_mean != is.na(impenv_mean)
    )

control_group <- 
  HANZE_and_ESS10 |> 
  filter(
    is.na(impenv_mean)
  )


#plot it on a map

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




