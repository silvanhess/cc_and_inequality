# Load Packages --------------------------------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(tmap)
library(ggplot2)


# Import data -----------------------------------------------------------------------------------------------------

# geographical data
# download from https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics
# NUTS 2010

sf_nuts <- read_sf("data/NUTS/NUTS_RG_20M_2010_3035.shp/NUTS_RG_20M_2010_3035.shp")

# survey data
# Download from https://ess.sikt.no/en/

ESS8 <- read_csv("data/ESS/ESS8e02_3/ESS8e02_3.csv")
ESS10 <- bind_rows(read_csv("data/ESS/ESS10/ESS10.csv"), read_csv("data/ESS/ESS10SC/ESS10SC.csv"))
ESS11 <- read_csv("data/ESS/ESS11/ESS11.csv")
ESS_all <- read_csv("data/ESS/ESSsubset/ESS8e02_3-ESS9e03_2-ESS10-ESS10SC-ESS11-subset.csv")

# flood data
# download from https://essd.copernicus.org/articles/16/5145/2024/essd-16-5145-2024-assets.html

# HANZE <- read_csv("data/HANZE/HANZE_databsae_1870-2020_Version_v.2.1.2/HANZE_events.csv")
floods_regions <- read_sf("data/HANZE/HANZE_databsae_1870-2020_Version_v.2.1.2/HANZE_floods_regions_2010/HANZE_floods_regions_2010.shp")


# Prepare Geographical Data ---------------------------------------------------------------------------------------

# prepare a dataset with geographical information of europe

nuts <- 
  sf_nuts |> 
  rename(region = NUTS_ID) |> 
  dplyr::select(region, CNTR_CODE, LEVL_CODE, NAME_LATN, geometry)

colnames(sf_nuts)
class(sf_nuts)


country_polygons <- nuts |> filter(LEVL_CODE == 0)

# Prepare survey data ---------------------------------------------------------------------------------------------

# variables
# impenv, importance to care for environment
# wrclmch, how worried about climate change
# ccnthum, CC caused naturally or by humans
# ccrdprs, personal responsability to reduce CC

# prepare a dataset with survey participants from wave eight and ten

ESS_prepared <- 
  bind_rows(ESS10, ESS8) |> 
  filter(essround == 8 | essround == 10) |> 
  filter(cntry %in% ESS8$cntry, # make sure that the country in wave 10 is also in wave 8
         cntry %in% ESS10$cntry, # make sure that the country in wave 8 is also in wave 10
         cntry %in% c("GB", "ES", "CH", "AT") # only include countries with sufficient data
   ) |> 
  mutate(
    date8 = as.Date(paste(inwdds, inwmms, inwyys, sep = "-"), format = "%d-%m-%Y"),
    date10 = inwds
  ) |>
  mutate(
    date = if_else(essround == 8, date8, date10) |> as.Date(),
    cntry = if_else(cntry == "GB", "UK", cntry)
  ) |>
  left_join(y = nuts, by = "region") |>
  dplyr::select(idno, essround, date, 
         cntry, region, LEVL_CODE, NAME_LATN, # geographical information
         wrclmch, # outcome variable
         lrscale, ccnthum, # control variables
         geometry) |>
  mutate(
    region = if_else(region == "99999", NA, region),
    wrclmch = if_else(wrclmch > 5, NA, wrclmch),
    lrscale = if_else(lrscale > 10, NA, lrscale),
    ccntum = if_else(ccnthum > 5, NA, ccnthum)
  ) |>
  drop_na(idno, essround, cntry, region) |>
  mutate(respondent_id = row_number()) |> 
  dplyr::select(respondent_id, essround, date, 
         cntry, region, LEVL_CODE, NAME_LATN, # geographical information
         wrclmch, # outcome variable
         lrscale, ccnthum, # control variables
         geometry) |> 
  st_as_sf()

# ESS_prepared |>
#   as_tibble() |>
#   group_by(cntry) |>
#   summarize(count = n()) |>
#   print(n=999)

# ESS_prepared |>
#   as_tibble() |> 
#   filter(essround == 10) |>
#   group_by(cntry) |>
#   summarise(min_date = min(date, na.rm = TRUE),
#             mean_date = mean(date, na.rm = TRUE),
#             max_date = max(date, na.rm = TRUE))



# Prepare Floods Data ---------------------------------------------------------------------------------------------

# prepare a dataset of floods from 2016 until 2019 in the areas where people were surveyed

floods_prepared <- 
  floods_regions |> 
  mutate(end_date = make_date(End_Y, End_M, End_D)) |> 
  filter(between(as.Date(end_date, format = "%Y-%m-%d"), 
                 as.Date("2017-06-18"), # date of the last survey in wave 8
                 as.Date("2021-05-05")), # date of the first survey in wave 10
         Code %in% ESS_prepared$cntry # only floods in countries that were surveyed
         ) |>
  rename(
    country = Code,
    flood_id = ID,
    regions = Region2010
  ) |> 
  dplyr::select(flood_id, end_date, country, regions)
  

# Visualize Survey Regions and Floods  ---------------------------------------------------------------------------------------------------

ESS_prepared_regions <- 
  ESS_prepared |> 
  group_by(region, geometry) |> 
  summarise(count = n())

map <- tm_basemap("OpenStreetMap") +
  tm_shape(floods_prepared$geometry) +
  tm_polygons(fill = "blue", fill_alpha = 0.5) +
  tm_borders(col = "black", lwd = 1) +
  tm_shape(ESS_prepared_regions$geometry) +
  tm_polygons(fill = "red", fill_alpha = 0.5) +
  tm_borders(col = "black", lwd = 1)

tmap_mode("view")
# print(map)



# Spatial Merge ---------------------------------------------------------------------------------------------------

sf_ESS_prepared_grouped <- 
  ESS_prepared |> 
  group_by(region, geometry) |> 
  summarise(count = n())

containment_check_old <-
  st_intersects(sf_ESS_prepared_grouped, floods_prepared)

# Shrink each polygon by 10 kilometers (to avoid overlaps to neighbouring regions)
sf_ESS_shrunk <- st_buffer(sf_ESS_prepared_grouped, dist = -10000)

containment_check <- 
  st_intersects(sf_ESS_shrunk, floods_prepared)

ESS_intersects <- 
  st_join(sf_ESS_shrunk, floods_prepared, join = st_intersects) |> 
  distinct(region, .keep_all = TRUE) |>  # to flag a region as flooded only one flood is necessary
  as_tibble()

ESS_prepared_with_flood_info <- 
  left_join(ESS_prepared, ESS_intersects, by = "region") |> 
  mutate(flood = if_else(!is.na(flood_id), 1, 0)) |> 
  dplyr::select(respondent_id, essround, date, cntry, region, LEVL_CODE, NAME_LATN, 
                wrclmch, lrscale, ccnthum, geometry.x, flood) |> 
  rename(geometry = geometry.x)



# Visualize Treatment Group & Control Group -----------------------------------------------------------------------

treatment_region_before <- 
  ESS_prepared_with_flood_info |> 
  filter(flood == 1,
         essround == 8) |> 
  group_by(region, geometry) |> 
  summarise(count = n())

treatment_region_after <- 
  ESS_prepared_with_flood_info |> 
  filter(flood == 1,
         essround == 10) |> 
  group_by(region, geometry) |> 
  summarise(count = n())

control_region_before <- 
  ESS_prepared_with_flood_info |> 
  filter(flood == 0,
         essround == 8) |> 
  group_by(region, geometry) |> 
  summarise(count = n())

control_region_after <- 
  ESS_prepared_with_flood_info |> 
  filter(flood == 0,
         essround == 10) |> 
  group_by(region, geometry) |> 
  summarise(count = n())

map <- tm_basemap("OpenStreetMap") +
  tm_shape(floods_prepared$geometry) +
  tm_polygons(fill = "blue", fill_alpha = 0.5) +
  tm_borders(col = "black", lwd = 1) +  
  tm_shape(treatment_region_before$geometry) +
  tm_polygons(fill = "red", fill_alpha = 0.5) +
  tm_borders(col = "black", lwd = 1) + 
  tm_shape(treatment_region_after$geometry) +
  tm_polygons(fill = "green", fill_alpha = 0.5) +
  tm_borders(col = "black", lwd = 1) +
  tm_shape(control_region_before$geometry) +
  tm_polygons(fill = "purple", fill_alpha = 0.5) +
  tm_borders(col = "black", lwd = 1) + 
  tm_shape(control_region_after$geometry) +
  tm_polygons(fill = "yellow", fill_alpha = 0.5) +
  tm_borders(col = "black", lwd = 1) +
  tm_shape(country_polygons$geometry) +
  tm_polygons(fill = "orange", fill_alpha = 0.5) +
  tm_borders(col = "black", lwd = 1)

tmap_mode("view")
# print(map)


# Code the Treatment Variable & Outcome Variable -------------------------------------------------------------------

data <- 
  ESS_prepared_with_flood_info |> 
  as_tibble() |> 
  mutate(time = if_else(essround == 8, 0, 1)) |> 
  rename(
    treatment_variable = flood,
    outcome_variable = wrclmch
    # control_variable = lrscale
  ) |> 
  mutate(
    outcome_variable = factor(outcome_variable, ordered = TRUE)
    # control_variable = factor(control_variable, ordered = TRUE)
  )

# Check Pre-Treatment Characteristics -----------------------------------------------------------------------------


pre_treatment_data <- data |>
  filter(time == 0)  # Keep only rows from the pre-treatment period

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

logit_model <- glm(outcome_variable ~ time + treatment_variable + interaction, data = data, family = "binomial")
summary(logit_model)

library(MASS)
ordered_logit_model <- polr(outcome_variable ~ time + treatment_variable + interaction, data = data, method = "logistic")
summary(ordered_logit_model)


# Robustness Checks -----------------------------------------------------------------------------------------------

# 1. Multinomial Logit Model

library(nnet)
mnl_model <- multinom(outcome_variable ~ time + treatment_variable + interaction, data = data)
summary(mnl_model)

# 2. Parallel Trends Assumption Check

# data$placebo_time <- ifelse(data$time < actual_treatment_time, 1, 0)
# placebo_model <- glm(outcome_variable ~ placebo_time * treatment_variable, 
#                      family = "binomial", data = data)
# summary(placebo_model)

# 3. Heterogeneous Treatment Effects

heterogeneity_model1 <- glm(outcome_variable ~ time + treatment_variable + interaction * lrscale, 
                           family = "binomial", data = data)
summary(heterogeneity_model1)

heterogeneity_model2 <- glm(outcome_variable ~ time + treatment_variable + interaction * ccnthum, 
                            family = "binomial", data = data)
summary(heterogeneity_model2)

# --> treatment effect does not vary by subgroup

# 4. Visualization of Treatment Effects

# Outcome distribution before vs. after treatment
ggplot(data, aes(x = outcome_variable, fill = as.factor(time))) +
  geom_bar(position = "dodge") +
  labs(title = "Outcome Distribution Before & After Treatment") +
  facet_wrap(facets = ~cntry)



