#Author: Nicola Smith
#Date: 16 July 2025
#Project: Climate change enhances the success of marine invasive species
#Aim: Create invasion hotspot maps

getwd()
setwd("/Users/nicolasmith/Desktop/Nic")

library(ggplot2)
library(dplyr)

install.packages("sf")
library(sf)
install.packages("rnaturalearth")
library(rnaturalearth)
install.packages("rnaturalearthdata")
library(rnaturalearthdata)
install.packages("viridis")
library(viridis)
install.packages("viridisLite")
library(viridisLite)

theme_set(theme_minimal())

# The new versions do not have the EEZ boundary, they are cropped to the PPOW (pelagic provinces of the world)
# meow_ppow <- read_sf("../DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1//01_Data//WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp")
# meow <- meow_ppow |> filter(TYPE == "MEOW")

# Old version from here:
# https://web.archive.org/web/20120415041842/http://conserveonline.org/workspaces/ecoregional.shapefile/MEOW/

meow <- read_sf("MEOW2/meow_ecos.shp")
coast_fill <- rnaturalearth::ne_countries()
coast <- rnaturalearth::ne_coastline()

dat0 <- readr::read_csv("Donar region_Total_no_species.csv")

dat <- left_join(dat0, meow, by = c("Native_range_marine_province" = "PROVINCE")) |>
  st_as_sf()

ggplot() +
  # geom_sf(data = meow) +
  geom_sf(data = dat, aes(fill = no._species)) +
  geom_sf(data = coast_fill, fill = "white", colour = "grey85") +
  # geom_sf(data = coast, colour = "grey85") +
  viridis::scale_fill_viridis(breaks = c(1, 3, 5, 7, 9)) + # I don't know what palette you used
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey95")) +
  coord_sf(expand = FALSE) +
  labs(fill = "No. of species")


meow <- read_sf("MEOW2/meow_ecos.shp")
coast_fill <- rnaturalearth::ne_countries()
coast <- rnaturalearth::ne_coastline()

dat1 <- readr::read_csv("Recipient_regions_Total_no._species")

dat <- left_join(dat0, meow, by = c("Invaded_range_Marine_province" = "PROVINCE")) |>
  st_as_sf()

ggplot() +
  # geom_sf(data = meow) +
  geom_sf(data = dat, aes(fill = no._species)) +
  geom_sf(data = coast_fill, fill = "white", colour = "grey85") +
  # geom_sf(data = coast, colour = "grey85") +
  viridis::scale_fill_viridis(breaks = c(1, 3, 5, 7, 9,11,13,15,17)) + # I don't know what palette you used
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey95")) +
  coord_sf(expand = FALSE) +
  labs(fill = "No. of species")
