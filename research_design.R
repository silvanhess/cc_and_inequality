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

# class(sf_nuts)
# sf_nuts$NUTS_ID

nuts <- 
  sf_nuts |> 
  as_tibble() |> 
  rename(region = NUTS_ID) |> 
  select(region,LEVL_CODE, NAME_LATN, geometry)

# nuts |> 
#   filter(str_starts(region, "IT")) |> 
#   view()

# Prepare dataset ---------------------------------------------------------------------------------------------

# prepare a dataset with survey participants from wave eight and ten

ESS <- bind_rows(ESS10, ESS8)

ESS_prepared <-
  left_join(x=ESS, y=nuts) |> 
  select(idno, essround, cntry, region, LEVL_CODE, NAME_LATN , agea, gndr,lrscale, impenv, geometry) |> 
  drop_na(idno, essround, cntry, region, agea, gndr,lrscale, impenv) |>
  filter(
    region != "99999",
    lrscale <= 10,
    impenv <= 6
  ) |> 
  mutate(
    gndr = factor(gndr, levels = c(1, 2), labels = c("Male", "Female"))
  )

# prepare a dataset of floods from 2016 until 2019 in the areas where people were surveyed

HANZE_events_prepared <- 
  HANZE_events |>
    rename(
      end_date = `End date`,
      country = `Country code`,
      region = `Regions affected (v2010)`
    ) |> 
    select(ID, end_date, country, region) |> 
    filter(
      between(as.Date(end_date, format = "%Y-%m-%d"), as.Date("2016-01-01"), as.Date("2019-12-31")) # maybe use year variable instead
    ) |> 
    separate_rows(region, sep = ";")


HANZE_events_prepared_geom <- left_join(x=HANZE_events_prepared, y = nuts) |>
  drop_na(LEVL_CODE)

level_2_countries <- c("AT", "BE", "CH", "ES", "FR", "GR", "IT", "NL", "NO", "PL", "PT")



# Spatial Merge from level 3 to level 2 ---------------------------------------------------------------------------------------------------

library(sf)
library(dplyr)

# Example level_2_countries vector
level_2_countries <- c("AT", "BE") # Replace with your actual Level 2 countries

# Filter HANZE_events_prepared_geom for Level 2 countries
HANZE_level2 <- HANZE_events_prepared_geom %>%
  filter(country %in% level_2_countries)

# Filter nuts dataset to include only Level 2 regions
nuts_level2 <- nuts %>%
  filter(LEVL_CODE == 2)

# Perform spatial join to associate each region in HANZE_level2 with its Level 2 region in nuts
HANZE_level2_aggregated <- HANZE_level2 %>%
  st_join(nuts_level2, join = st_within, left = FALSE) %>%
  group_by(region = nuts_level2$region) %>%
  summarize(
    ID = first(ID), # Adjust aggregation for ID as needed
    end_date = max(end_date), # Adjust as per your aggregation rule
    NAME_LATN = first(NAME_LATN), # Adjust as per your aggregation rule
    geometry = st_union(geometry)
  )

# Combine with non-Level 2 data
HANZE_other <- HANZE_events_prepared_geom %>%
  filter(!country %in% level_2_countries)

# Combine aggregated Level 2 data and other data
HANZE_combined <- bind_rows(HANZE_level2_aggregated, HANZE_other)

# Save or inspect the result
st_write(HANZE_combined, "HANZE_aggregated.shp")



# Create the Treatment and Control Groups -------------------------------------------------------------------------




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



# Plot the events over time ---------------------------------------------------------------------------------------

# Plot histogram
ggplot(HANZE_events_prepared, aes(x = end_date)) +
  geom_histogram(binwidth = 30, fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of End Dates",
    x = "End Date",
    y = "Count"
  ) +
  theme_minimal()


# Create map ------------------------------------------------------------------------------------------------------

ESS_prepared_geom_3 <- 
  ESS_prepared |> 
  group_by(region, cntry, LEVL_CODE, geometry) |> 
  summarise(count = n()) |> 
  filter(LEVL_CODE == 3)

ESS_prepared_geom_3 |> 
  group_by(cntry) |> 
  summarise(count = n())

ESS_prepared_geom_2 <- 
  ESS_prepared |> 
  group_by(region, cntry, LEVL_CODE, geometry) |> 
  summarise(count = n()) |> 
  filter(LEVL_CODE == 2)

ESS_prepared_geom_2 |> 
  group_by(cntry) |> 
  summarise(count = n()) |> 
  select(cntry) |> 
  as.character()

ESS_prepared_geom_1 <- 
  ESS_prepared |> 
  group_by(region, cntry, LEVL_CODE, geometry) |> 
  summarise(count = n()) |> 
  filter(LEVL_CODE == 1)

ESS_prepared_geom_1 |> 
  group_by(cntry) |> 
  summarise(count = n())

HANZE_geom <- 
  left_join(x = HANZE_events_prepared, y = nuts) |> 
  drop_na(LEVL_CODE) |> 
  group_by(region, LEVL_CODE, geometry) |> 
  summarise(count = n())

treatment_group_geom <- treatment_group |> 
  group_by(region, LEVL_CODE, geometry) |> 
  summarise(count = n())

map <- tm_basemap("OpenStreetMap") +
  tm_shape(ESS_prepared_geom_3$geometry) +
  tm_polygons(col = "red") +
  tm_shape(ESS_prepared_geom_2$geometry) +
  tm_polygons(col = "green") +
  tm_shape(ESS_prepared_geom_1$geometry) +
  tm_polygons(col = "purple") +
  tm_shape(HANZE_geom$geometry) +
  tm_polygons(col = "blue") +
  tm_shape(treatment_group_geom$geometry) +
  tm_polygons(col = "yellow")
tmap_mode("view")
tmap_options(check.and.fix = TRUE)
print(map)





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


