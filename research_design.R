# Load Libraries --------------------------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(haven)
library(labelled)
library(countrycode)
library(tmap)
library(sf) # für shape files


# Import data -----------------------------------------------------------------------------------------------------

# deeznuts

sf_nuts <- read_sf("data/NUTS/NUTS_RG_20M_2010_3035.shp/NUTS_RG_20M_2010_3035.shp")

# survey data

ESS10 <- read_csv("data/ESS/ESS10/ESS10.csv")
ESS8 <- read_csv("data/ESS/ESS8e02_3/ESS8e02_3.csv")

# flood data

HANZE <- read_csv("data/HANZE/v2.1/HANZE_events.csv")


# Prepare dataset ---------------------------------------------------------------------------------------------

# prepare a dataset with geographical information of europe

nuts <- 
  sf_nuts |> 
  as_tibble() |> 
  rename(region = NUTS_ID) |> 
  select(region,LEVL_CODE, NAME_LATN, geometry)

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
    gndr = factor(gndr, levels = c(1, 2), labels = c("Male", "Female"))
  )

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
  filter(LEVL_CODE == 2) |> 
  st_as_sf()

# Perform spatial join to check if the events are within the Level 2 regions
containment_check <- st_within(HANZE_level2_old, nuts_level2)

# Join the small polygons with the large ones based on containment
joined <- st_join(HANZE_level2_old, nuts_level2, join = st_within)

# select the columns of the level 2 regions
HANZE_level2_new <- 
  joined |> 
    select(ID, end_date, country, region.y, LEVL_CODE.y, NAME_LATN.y, geometry) |>
    rename(region = region.y, LEVL_CODE = LEVL_CODE.y, NAME_LATN = NAME_LATN.y)
  

# Spatial Merge from level 3 to level 1 --------------------------------------------------------------------------------

level_1_countries <- c("DE", "GB", "ME") # countries that only have level 1 regions in the ESS dataset

# Filter HANZE_events_prepared_geom for countries that have Level 2 regions
HANZE_level1_old <- 
  HANZE_prepared %>%
  filter(country %in% level_1_countries) |> 
  st_as_sf()

# Filter nuts dataset to include only Level 1 regions
nuts_level1 <- 
  nuts %>%
  filter(LEVL_CODE == 1) |> 
  st_as_sf()

# Perform spatial join to check if the events are within the Level 1 regions
containment_check <- st_within(HANZE_level1_old, nuts_level1)

# Join the small polygons with the large ones based on containment
joined <- st_join(HANZE_level1_old, nuts_level1, join = st_within)

# select the columns of the level 1 regions
HANZE_level1_new <- 
  joined |> 
  select(ID, end_date, country, region.y, LEVL_CODE.y, NAME_LATN.y, geometry) |>
  rename(region = region.y, LEVL_CODE = LEVL_CODE.y, NAME_LATN = NAME_LATN.y)


# Put the three Datasets together ---------------------------------------------------------------------------------

# filter countries that are neither in level_1_countries nor in level_2_countries
HANZE_level3 <- 
  HANZE_prepared |> 
  filter(!(country %in% level_1_countries) & !(country %in% level_2_countries))

HANZE_combined <- 
  bind_rows(HANZE_level3, HANZE_level2_new, HANZE_level1_new) |> 
  arrange(ID) |> 
  distinct(across(-geometry), .keep_all = TRUE) # remove duplicates

HANZE_combined

# Create the Treatment and Control Groups -------------------------------------------------------------------------

# identify the surveyed people that experienced floods
# put them in the treatment group

treatment_group <-     
  left_join(x = HANZE_combined, 
            y = ESS_prepared
  ) |> 
  group_by(across(-region)) |> 
  summarize(
    region = paste(unique(region), collapse = ";"),
    .groups = "drop")


treatment_group_before <- treatment_group |> filter(essround == 8)
treatment_group_after <- treatment_group |> filter(essround == 10) 


# identify the surveyed people that haven't experienced floods
# put them in the control group


control_group <- ESS_prepared %>%
  filter(!idno %in% treatment_group$idno)

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

treatment_group$ID |> n_distinct() # 96 Floods in the treatment group
treatment_group$region |> n_distinct() # 136 affected regions in the treatment group

treatment_group |> 
  group_by(ID) |> 
  summarize(count = n()) |> 
  view() # many Floods only have 1 Person -> drop them


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
    LEVL_CODE == 2,
    cntry %in% level_2_countries)

ESS_prepared_geom_1 <- 
  ESS_prepared |> 
  group_by(region, cntry, LEVL_CODE, geometry) |> 
  summarise(count = n()) |> 
  filter(cntry %in% level_1_countries)


HANZE_geom <- 
  left_join(x = HANZE_combined, y = nuts) |> 
  drop_na(LEVL_CODE) |> 
  group_by(region, LEVL_CODE, geometry) |> 
  summarise(count = n())

treatment_group_geom <- treatment_group |> 
  group_by(region, LEVL_CODE, geometry) |> 
  summarise(count = n())

control_group_geom <- control_group |> 
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
  tm_polygons(col = "blue")
  # tm_shape(treatment_group_geom$geometry) +
  # tm_polygons(col = "yellow") +
  # tm_shape(control_group$geometry) +
  # tm_polygons(col = "black")
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


