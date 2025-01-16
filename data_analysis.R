# Load Libraries --------------------------------------------------------------------------------------------------

library(tidyverse)
# library(readxl)
# library(haven)
# library(labelled)
# library(countrycode)
library(sf) # für shape files
library(tmap)



# Import data -----------------------------------------------------------------------------------------------------

# geographical data
# download from https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics
# NUTS 2010

sf_nuts <- read_sf("data/NUTS/NUTS_RG_20M_2010_3035.shp/NUTS_RG_20M_2010_3035.shp")

# Alternative: use the files from the HANZE database
# read_sf("data/HANZE/HANZE_databsae_1870-2020_Version_v.2.1.2/HANZE_floods_regions_2010/HANZE_floods_regions_2010.shp")

# sf_nuts |> 
#   group_by(CNTR_CODE) |>
#   summarize(count = n()) |>
#   view()

# survey data
# Download from https://ess.sikt.no/en/

ESS10 <- read_csv("data/ESS/ESS10/ESS10.csv")
ESS10SC <- read_csv("data/ESS/ESS10SC/ESS10SC.csv")
ESS10 <- bind_rows(ESS10, ESS10SC)
ESS8 <- read_csv("data/ESS/ESS8e02_3/ESS8e02_3.csv")

# ESS10 |> filter(cntry == "IT") |> select(regunit)
# ESS8 |> filter(cntry == "IT") |> select(regunit)
# whiy does italy have level 1 regions in wave 10??

# flood data
# download from https://essd.copernicus.org/articles/16/5145/2024/essd-16-5145-2024-assets.html

HANZE <- read_csv("data/HANZE/HANZE_databsae_1870-2020_Version_v.2.1.2/HANZE_events.csv")


# Prepare dataset ---------------------------------------------------------------------------------------------

# prepare a dataset with geographical information of europe

nuts <- 
  sf_nuts |> 
  as_tibble() |> 
  rename(region = NUTS_ID) |> 
  select(region, CNTR_CODE, LEVL_CODE, NAME_LATN, geometry)



# prepare a dataset with survey participants from wave eight and ten

ESS_prepared <- 
  bind_rows(ESS10, ESS8) |> 
  filter(region %in% ESS8$region, # make sure that the region in wave 10 is also in wave 8
         region %in% ESS10$region) |> # make sure that the region in wave 8 is also in wave 10
  mutate(
    date8 = as.Date(paste(inwdds, inwmms, inwyys, sep = "-"), format = "%d-%m-%Y"),
    date10 = inwds
  ) |> 
  mutate(
    date = if_else(essround == 8, date8, date10) |> as.Date()
  ) |> 
  left_join(y=nuts) |> 
  select(idno, essround, date, cntry, region, LEVL_CODE, NAME_LATN, lrscale, wrclmch, geometry) |> 
  drop_na(idno, essround, cntry, region,  lrscale, wrclmch) |>
  filter(
    region != "99999",
    lrscale <= 10,
    wrclmch <= 5
  ) |> 
  # mutate(
  #   cntry = if_else(cntry == "GB", "UK", cntry) # change GB to UK in cntry
  # ) |> 
  filter(LEVL_CODE != 1) |>  # exclude level 1 regions
  mutate(respondent_id = seq_len(n())) |> 
  select(respondent_id, essround, date, cntry, region, LEVL_CODE, NAME_LATN, agea, gndr, lrscale, wrclmch, geometry)

# ESS_prepared |>
#   group_by(region, essround) |>
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
      region = `Regions affected (v2010)`,
      flood_id = ID
    ) |> 
    select(flood_id, end_date, country, region) |> 
    filter(
      between(as.Date(end_date, format = "%Y-%m-%d"), as.Date("2017-11-01"), as.Date("2021-02-01"))
    ) |> 
    separate_rows(region, sep = ";") |> 
    left_join(y = nuts) |> 
    drop_na(LEVL_CODE) |> 
    arrange(flood_id) |> 
    distinct(across(-geometry), .keep_all = TRUE)

# HANZE_prepared |> 
#   group_by(ID) |> 
#   summarise(count = n())

# Spatial Merge from level 3 to level 2 ---------------------------------------------------------------------------------

# Identify Countries that have level 2 regions
ESS_prepared |> 
  filter(LEVL_CODE == 2) |> 
  group_by(cntry) |> 
  summarise(count = n()) 
level_2_countries <- c("AT", "BE", "CH", "ES", "FR", "NL", "NO", "PL", "PT") 

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
  select(flood_id, end_date, country, region.y, LEVL_CODE.y, NAME_LATN.y) |>
  as_tibble() |>  
  select(-geometry) |> 
  rename(region = region.y, LEVL_CODE = LEVL_CODE.y, NAME_LATN = NAME_LATN.y) |> 
  distinct(flood_id, region, .keep_all = TRUE) |> # delete duplicates
  left_join(y = nuts) |> # add the polygons of level 2 regions
  select(-CNTR_CODE)


  
# visualize on map
# map <- tm_basemap("OpenStreetMap") +
#   tm_shape(HANZE_level2_new$geometry) +
#   tm_polygons(col = "black")
# tmap_mode("view")
# print(map)


# Spatial Merge from level 3 to level 1 --------------------------------------------------------------------------------

# # Identify Countries that have level 1 regions
# ESS_prepared |> 
#   filter(LEVL_CODE == 1) |> 
#   group_by(cntry) |> 
#   summarise(count = n()) 
# 
# level_1_countries <- c("DE", "UK")
# 
# # Filter HANZE_events_prepared_geom for countries that have Level 2 regions
# HANZE_level1_old <- 
#   HANZE_prepared %>%
#   filter(country %in% level_1_countries) |> 
#   st_as_sf()
# 
# # Filter nuts dataset to include only Level 1 regions
# nuts_level1 <- 
#   nuts %>%
#   filter(LEVL_CODE == 1,
#          CNTR_CODE %in% level_1_countries) |> 
#   st_as_sf()
# 
# # Perform spatial join to check if the events are within the Level 1 regions
# containment_check <- st_intersects(HANZE_level1_old, nuts_level1)
# 
# # Join the small polygons with the large ones based on containment
# level1_joined <- st_join(HANZE_level1_old, nuts_level1, join = st_intersects)
# 
# # select the columns of the level 1 regions
# HANZE_level1_new <- 
#   level1_joined |> 
#   select(ID, end_date, country, region.y, LEVL_CODE.y, NAME_LATN.y) |>
#   as_tibble() |>  
#   select(-geometry) |> # remove the polygons of level 3 regions
#   rename(region = region.y, LEVL_CODE = LEVL_CODE.y, NAME_LATN = NAME_LATN.y) |> 
#   distinct(ID, region, .keep_all = TRUE) |> # delete duplicates
#   left_join(y = nuts) |> # add the polygons of level 1 regions
#   select(-CNTR_CODE)
# 
# # visualize on map
# map <- tm_basemap("OpenStreetMap") +
#   tm_shape(HANZE_level2_new$geometry) +
#   tm_polygons(col = "black") +
#   tm_shape(HANZE_level1_new$geometry) +
#   tm_polygons(col = "red")
# # tmap_mode("view")
# # print(map)


# Put the Datasets together ---------------------------------------------------------------------------------

HANZE_level3 <- 
  HANZE_prepared |> 
  filter(!(country %in% level_2_countries))

HANZE_combined <- 
  bind_rows(HANZE_level3, HANZE_level2_new) |> 
  mutate(
    LEVL_CODE = if_else(is.na(LEVL_CODE), 2, LEVL_CODE)
  ) |> 
  select(flood_id, end_date, country, region, LEVL_CODE, NAME_LATN, geometry)

# HANZE_combined |> 
#   group_by(country) |> 
#   summarise(count = n())

# visualize on map
# map <- tm_basemap("OpenStreetMap") +
#   tm_shape(HANZE_combined$geometry) +
#   tm_polygons(col = "black")
#tmap_mode("view")
#print(map)


# Create the Treatment and Control Groups -------------------------------------------------------------------------

# identify the surveyed people that experienced floods
# put them in the treatment group

treatment_group <-     
  inner_join(x = ESS_prepared, 
            y = HANZE_combined,
            relationship = "many-to-many"
  ) |> 
  group_by(across(-region)) |> 
  summarize(
    region = paste(unique(region), collapse = ";"),
    .groups = "drop") |> 
  select(
    respondent_id, essround, date, cntry, region, LEVL_CODE, NAME_LATN, agea, gndr, lrscale, wrclmch, flood_id, end_date, geometry
  ) |>  # arrange the variables manually
  distinct(respondent_id, .keep_all = TRUE) # delete duplicates

# treatment_group |> 
#   group_by(region) |> 
#   summarize(count = n()) |> 
#   view() # some regions only have 1 Person -> drop them
# 
# region_selector <- 
#   treatment_group |> 
#     group_by(region) |> 
#     summarize(count = n()) |> 
#     filter(count >= 30) |>
#     select(region) |> 
#     as.list()
# 
# 
# # create dataframe with regions that have at least 30 obvervations
# treatment_group_filtered <- 
#   treatment_group |> 
#     filter(region %in% region_selector$region)
# 
# 
# treatment_group_filtered |> 
#   group_by(region) |> 
#   summarize(count = n()) |> 
#   view() # now we have regions with at least 30 observations



# split into before and after the floodings
treatment_group_before <- treatment_group |> filter(essround == 8)
treatment_group_after <- treatment_group |> filter(essround == 10) 


# identify the surveyed people that haven't experienced floods
# put them in the control group

control_group <- ESS_prepared |> 
  anti_join(treatment_group, by = "respondent_id") |> 
  filter(cntry %in% treatment_group$cntry)

# split into before and after the floodings
control_group_before <- control_group |> filter(essround == 8)
control_group_after <- control_group |> filter(essround == 10) 
  

# visualize on a map
# treatment_group_geom <- treatment_group |>
#   group_by(region, LEVL_CODE, geometry) |>
#   summarise(count = n())
# control_group_geom <- control_group |>
#   group_by(region, LEVL_CODE, geometry) |>
#   summarise(count = n())
# map <- tm_basemap("OpenStreetMap") +
#   tm_shape(treatment_group_geom$geometry) +
#   tm_polygons(col = "black") +
#   tm_tiles() +
#   tm_shape(control_group_geom$geometry) +
#   tm_polygons(col = "red")
# tmap_mode("view")
# print(map)



# matrix with number of survey respondents
survey_respondents_per_group <- 
  matrix(
    data = c(
      nrow(treatment_group_before),
      nrow(treatment_group_after),
      nrow(treatment_group),
      nrow(control_group_before),
      nrow(control_group_after),
      nrow(control_group),
      nrow(treatment_group_before) + nrow(control_group_before),
      nrow(treatment_group_after) + nrow(control_group_after),
      nrow(treatment_group) + nrow(control_group)
    ),
    nrow = 3,
    ncol = 3,
    dimnames = list(c("Treatment Group", "Control Group", "Total"), c("Before", "After", "Total")),
    byrow = TRUE
  )

# Descriptive Analysis ---------------------------------------------------------------------------------------

# Distribution of Floods over time
df_hist_HANZE_combined <- 
  HANZE_combined |> 
  distinct(ID, .keep_all = TRUE)
ggplot(df_hist_HANZE_combined, aes(x = end_date)) +
  geom_histogram(binwidth = 30, fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of End Dates",
    x = "End Date",
    y = "Count"
  ) +
  theme_minimal()

# Distribution of Countries
treatment_group |> 
  group_by(cntry) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = reorder(cntry, count), y = count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of Countries in the Treatment Group",
    x = "Country",
    y = "Count"
  ) +
  theme_minimal()

control_group |> 
  group_by(cntry) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = reorder(cntry, count), y = count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(
    title = "Distribution of Countries in the Control Group",
    x = "Country",
    y = "Count"
  ) +
  theme_minimal()


# Code the Treatment Variable & Outcome Variable -------------------------------------------------------------------

row_number(ESS_prepared)

data <- 
  bind_rows(treatment_group, control_group) |> 
    mutate(
      Flood = if_else(respondent_id %in% treatment_group$respondent_id, 1, 0)
    ) |> 
  mutate(time = if_else(essround == 8, 0, 1)) |> 
  rename(
    treatment_variable = Flood,
    outcome_variable = wrclmch,
    control_variable = lrscale
  ) |> 
  mutate(
    outcome_variable = factor(outcome_variable, ordered = TRUE),
    control_variable = factor(control_variable, ordered = TRUE))


# data$outcome_variable
# data$control_variable

# Check Pre-Treatment Characteristics -----------------------------------------------------------------------------


pre_treatment_data <- data %>%
  filter(time == 0)  # Keep only rows from the pre-treatment period

pre_treatment_data %>%
  group_by(treatment_variable) %>%
  summarize(
    mean_outcome = mean(outcome_variable, na.rm = TRUE),
    median_outcome = median(outcome_variable, na.rm = TRUE),
    sd_outcome = sd(outcome_variable, na.rm = TRUE),
    count = n()
  )
wilcox_test <- wilcox.test(
  outcome_variable ~ treatment_variable,
  data = pre_treatment_data
)
print(wilcox_test)
library(ggplot2)
ggplot(pre_treatment_data, aes(x = factor(treatment_variable), y = outcome_variable)) +
  geom_boxplot() +
  labs(x = "Treatment Group", y = "Outcome Variable", title = "Pre-Treatment Distribution")
ggplot(pre_treatment_data, aes(x = outcome_variable, fill = factor(treatment_variable))) +
  geom_density(alpha = 0.5) +
  labs(x = "Outcome Variable", fill = "Treatment Group", title = "Pre-Treatment Density Plot")


# Calculate the DiD Interaction Term ------------------------------------------------------------------------------

data <- data %>%
  mutate(interaction = time * treatment_variable)


# Run the Regression ----------------------------------------------------------------------------------------------

model <- glm(outcome_variable ~ time + treatment_variable + interaction, data = data, family = "binomial")
summary(model)

# library(MASS)
# model <- polr(outcome_variable ~ time + treatment_variable + interaction, data = data, method = "logistic")
# summary(model)


