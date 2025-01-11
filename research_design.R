# Load Libraries --------------------------------------------------------------------------------------------------

library(tidyverse)
# library(readxl)
# library(haven)
# library(labelled)
# library(countrycode)
library(tmap)
library(sf) # für shape files


# Import data -----------------------------------------------------------------------------------------------------

# geographical data
# download from https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics
# NUTS 2010

sf_nuts <- read_sf("data/NUTS/NUTS_RG_20M_2010_3035.shp/NUTS_RG_20M_2010_3035.shp")

# sf_nuts |> 
#   group_by(CNTR_CODE) |>
#   summarize(count = n()) |>
#   view()

# survey data
# Download from https://ess.sikt.no/en/

ESS10 <- read_csv("data/ESS/ESS10/ESS10.csv")
ESS8 <- read_csv("data/ESS/ESS8e02_3/ESS8e02_3.csv")

# flood data
# download from https://essd.copernicus.org/articles/16/5145/2024/essd-16-5145-2024-assets.html

HANZE <- read_csv("data/HANZE/v2.1/HANZE_events.csv")


# Prepare dataset ---------------------------------------------------------------------------------------------

# prepare a dataset with geographical information of europe

nuts <- 
  sf_nuts |> 
  as_tibble() |> 
  rename(region = NUTS_ID) |> 
  select(region, CNTR_CODE, LEVL_CODE, NAME_LATN, geometry)



# prepare a dataset with survey participants from wave eight and ten

ESS <- 
  bind_rows(ESS10, ESS8) |> 
  drop_na(cntry) |> 
  filter(cntry %in% ESS10$cntry) |>   # make sure that the country in wave 8 is also in wave 10
  mutate(
    date8 = as.Date(paste(inwdds, inwmms, inwyys, sep = "-"), format = "%d-%m-%Y"),
    date10 = inwds
  ) |> 
  mutate(
    date = if_else(essround == 8, date8, date10) |> as.Date()
  )

ESS_prepared <-
  left_join(x=ESS, y=nuts) |> 
  select(idno, essround, date, cntry, region, LEVL_CODE, NAME_LATN , agea, gndr,lrscale, impenv, geometry) |> 
  drop_na(idno, essround, cntry, region, agea, gndr,lrscale, impenv) |>
  filter(
    region != "99999",
    lrscale <= 10,
    impenv <= 6
  ) |> 
  mutate(
    gndr = factor(gndr, levels = c(1, 2), labels = c("Male", "Female")),
    cntry = if_else(cntry == "GB", "UK", cntry) # change GB to UK in cntry
  )

# ESS_prepared |> 
#   group_by(cntry) |> 
#   summarize(count = n()) |> 
#   view()

# ESS_prepared |> 
#   filter(essround == 8) |>
#   group_by(cntry) |>
#   summarise(min_date = min(date, na.rm = TRUE),
#             mean_date = mean(date, na.rm = TRUE),
#             max_date = max(date, na.rm = TRUE)) |> 
#   view()
# 
# ESS_prepared |> 
#   filter(essround == 10) |>
#   group_by(cntry) |>
#   summarise(min_date = min(date, na.rm = TRUE),
#             mean_date = mean(date, na.rm = TRUE),
#             max_date = max(date, na.rm = TRUE)) |> 
#   view()



# prepare a dataset of floods from 2016 until 2019 in the areas where people were surveyed

HANZE_prepared <- 
  HANZE |>
    rename(
      end_date = `End date`,
      country = `Country code`,
      region = `Regions affected (v2010)`
    ) |> 
    select(ID, end_date, country, region) |> 
    filter(
      between(as.Date(end_date, format = "%Y-%m-%d"), as.Date("2017-11-01"), as.Date("2021-02-01"))
    ) |> 
    separate_rows(region, sep = ";") |> 
    left_join(y = nuts) |> 
    drop_na(LEVL_CODE) |> 
    arrange(ID) |> 
    distinct(across(-geometry), .keep_all = TRUE)

# HANZE_prepared |> 
#   group_by(ID) |> 
#   summarise(count = n())

# Spatial Merge from level 3 to level 2 ---------------------------------------------------------------------------------

# countries that only have level 2 regions in the ESS dataset
level_2_countries <- c("AT", "BE", "CH", "ES", "FR", "GR", "IT", "NL", "NO", "PL", "PT") 

# Filter HANZE_events_prepared_geom for countries that have Level 2 regions
HANZE_level2_old <- 
  HANZE_prepared %>%
  filter(country %in% level_2_countries) |> 
  st_as_sf()

# Filter nuts dataset to include only Level 2 regions
nuts_level2 <- 
  nuts %>%
  filter(LEVL_CODE == 2,
         CNTR_CODE %in% level_2_countries) |> 
  st_as_sf()

# Perform spatial join to check if the events are within the Level 2 regions
containment_check <- st_intersects(HANZE_level2_old, nuts_level2)
# containment_check <- st_within(HANZE_level2_old, nuts_level2)

# Join the small polygons with the large ones based on containment
level2_joined <- st_join(HANZE_level2_old, nuts_level2, join = st_intersects)

# select the columns of the level 2 regions
HANZE_level2_new <- 
  level2_joined |> 
  select(ID, end_date, country, region.y, LEVL_CODE.y, NAME_LATN.y) |>
  as_tibble() |>  
  select(-geometry) |> 
  rename(region = region.y, LEVL_CODE = LEVL_CODE.y, NAME_LATN = NAME_LATN.y) |> 
  distinct(ID, region, .keep_all = TRUE) |> # delete duplicates
  left_join(y = nuts) |> # add the polygons of level 2 regions
  select(-CNTR_CODE)


  
# visualize on map
map <- tm_basemap("OpenStreetMap") +
  tm_shape(HANZE_level2_new$geometry) +
  tm_polygons(col = "black")
# tmap_mode("view")
# print(map)


# Spatial Merge from level 3 to level 1 --------------------------------------------------------------------------------

level_1_countries <- c("DE", "UK", "ME") # countries that only have level 1 regions in the ESS dataset

# Filter HANZE_events_prepared_geom for countries that have Level 2 regions
HANZE_level1_old <- 
  HANZE_prepared %>%
  filter(country %in% level_1_countries) |> 
  st_as_sf()

# Filter nuts dataset to include only Level 1 regions
nuts_level1 <- 
  nuts %>%
  filter(LEVL_CODE == 1,
         CNTR_CODE %in% level_1_countries) |> 
  st_as_sf()

# Perform spatial join to check if the events are within the Level 1 regions
containment_check <- st_intersects(HANZE_level1_old, nuts_level1)

# Join the small polygons with the large ones based on containment
level1_joined <- st_join(HANZE_level1_old, nuts_level1, join = st_intersects)

# select the columns of the level 1 regions
HANZE_level1_new <- 
  level1_joined |> 
  select(ID, end_date, country, region.y, LEVL_CODE.y, NAME_LATN.y) |>
  as_tibble() |>  
  select(-geometry) |> # remove the polygons of level 3 regions
  rename(region = region.y, LEVL_CODE = LEVL_CODE.y, NAME_LATN = NAME_LATN.y) |> 
  distinct(ID, region, .keep_all = TRUE) |> # delete duplicates
  left_join(y = nuts) |> # add the polygons of level 1 regions
  select(-CNTR_CODE)

# visualize on map
map <- tm_basemap("OpenStreetMap") +
  tm_shape(HANZE_level2_new$geometry) +
  tm_polygons(col = "black") +
  tm_shape(HANZE_level1_new$geometry) +
  tm_polygons(col = "red")
# tmap_mode("view")
# print(map)


# Put the three Datasets together ---------------------------------------------------------------------------------

# filter countries that are neither in level_1_countries nor in level_2_countries
HANZE_level3 <- 
  HANZE_prepared |> 
  filter(!(country %in% level_1_countries) & !(country %in% level_2_countries))

HANZE_combined <- 
  bind_rows(HANZE_level3, HANZE_level2_new, HANZE_level1_new) |> 
  mutate(
    LEVL_CODE = if_else(is.na(LEVL_CODE), 2, LEVL_CODE)
  )

# visualize on map
map <- tm_basemap("OpenStreetMap") +
  tm_shape(HANZE_combined$geometry) +
  tm_polygons(col = "black")
#tmap_mode("view")
#print(map)


# Create the Treatment and Control Groups -------------------------------------------------------------------------

# identify the surveyed people that experienced floods
# put them in the treatment group

treatment_group <-     
  left_join(x = HANZE_combined, 
            y = ESS_prepared,
            relationship = "many-to-many"
  ) |> 
  group_by(across(-region)) |> 
  summarize(
    region = paste(unique(region), collapse = ";"),
    .groups = "drop") |> 
  select(
    idno, essround, date, cntry, region, LEVL_CODE, NAME_LATN, agea, gndr, lrscale, impenv, ID, end_date, geometry
  ) # arrange the variables manually


treatment_group |> 
  group_by(region) |> 
  summarize(count = n()) |> 
  view() # some regions only have 1 Person -> drop them

region_selector <- 
  treatment_group |> 
    group_by(region) |> 
    summarize(count = n()) |> 
    filter(count >= 30) |>
    select(region) |> 
    as.list()


# create dataframe with regions that have at least 30 obvervations
treatment_group_filtered <- 
  treatment_group |> 
    filter(region %in% region_selector$region)


treatment_group_filtered |> 
  group_by(region) |> 
  summarize(count = n()) |> 
  view() # now we have regions with at least 30 observations


# visualize on a map
treatment_group_geom <- treatment_group_filtered |>
  group_by(region, LEVL_CODE, geometry) |>
  summarise(count = n())
map <- tm_basemap("OpenStreetMap") +
  tm_shape(treatment_group_geom$geometry) +
  tm_polygons(col = "black")
tmap_mode("view")
print(map)

treatment_group_before <- treatment_group_filtered |> filter(essround == 8)
treatment_group_after <- treatment_group_filtered |> filter(essround == 10) 


# identify the surveyed people that haven't experienced floods
# put them in the control group


control_group <- ESS_prepared %>%
  filter(!idno %in% treatment_group_filtered$idno)

control_group_before <- control_group |> filter(essround == 8)
control_group_after <- control_group |> filter(essround == 10) 


survey_respondents_per_group <- 
  matrix(
    data = c(
      nrow(treatment_group_before),
      nrow(treatment_group_after),
      nrow(treatment_group),
      nrow(control_group_before),
      nrow(control_group_after),
      nrow(control_group)
    ),
    nrow = 2,
    ncol = 3,
    dimnames = list(c("Treatment Group", "Control Group"), c("Before", "After", "Total")),
    byrow = TRUE
  )

# Plot the events over time ---------------------------------------------------------------------------------------

df_hist_HANZE_combined <- 
  HANZE_combined |> 
  distinct(ID, .keep_all = TRUE)

# Plot histogram
ggplot(df_hist_HANZE_combined, aes(x = end_date)) +
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

ESS_prepared_geom_2 <- 
  ESS_prepared |> 
  group_by(region, cntry, LEVL_CODE, geometry) |> 
  summarise(count = n()) |> 
  filter(
    LEVL_CODE == 2)

ESS_prepared_geom_1 <- 
  ESS_prepared |> 
  group_by(region, cntry, LEVL_CODE, geometry) |> 
  summarise(count = n()) |> 
  filter(
    LEVL_CODE == 1)

ESS_geom <- 
  ESS_prepared |> 
  group_by(region, LEVL_CODE, geometry) |> 
  summarize(number_of_participants_per_region = n(), .groups = "drop") |> 
  filter(
    LEVL_CODE >= 1)

HANZE_geom <- 
  left_join(x = HANZE_combined, y = nuts) |> 
  drop_na(LEVL_CODE) |> 
  group_by(region, LEVL_CODE, geometry) |> 
  summarise(count = n())

treatment_group_geom <- treatment_group_filtered |>
  group_by(region, LEVL_CODE, geometry) |>
  summarise(count = n())

control_group_geom <- control_group |>
  group_by(region, LEVL_CODE, geometry) |>
  summarise(count = n()) |> 
  drop_na(LEVL_CODE)

map <- tm_basemap("OpenStreetMap") +
  # tm_shape(ESS_prepared_geom_1$geometry) +
  # tm_polygons(col = "white") +
  # tm_shape(ESS_prepared_geom_2$geometry) +
  # tm_polygons(col = "white") +
  # tm_shape(ESS_prepared_geom_3$geometry) +
  # tm_polygons(col = "white") +
  # tm_shape(ESS_geom$geometry) +
  # tm_polygons(col = "red") +
  tm_shape(HANZE_geom$geometry) +
  tm_polygons(col = "blue") +
  tm_shape(treatment_group_geom$geometry) +
  tm_polygons(col = "yellow") +
  tm_shape(control_group_geom$geometry) +
  tm_polygons(col = "black")
tmap_mode("view")
tmap_options(check.and.fix = TRUE)
print(map)
