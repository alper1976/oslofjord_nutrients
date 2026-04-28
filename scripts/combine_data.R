library(tidyverse)
library(sf)
library(readxl)
library(lubridate)
library(dbscan)
library(ggspatial)


setwd("/Users/alexaei/Jottacloud/projects/CEDAP/data")

# --- 1. Load GBIF Archive (Improved Join Logic) ---
events <- read_tsv("event.txt")
measurements <- read_tsv("measurementorfacts.txt")

# --- Audit: See every unique parameter name in the raw GBIF data ---
raw_parameter_names <- measurements %>%
  distinct(measurementType) %>%
  arrange(measurementType)

# Print to console
print(raw_parameter_names, n = Inf) 

# --- Audit: See names paired with their units ---
# This helps identify if 'Phosphorus' is sometimes in 'µM' or 'mg/L'
type_unit_pairs <- measurements %>%
  distinct(measurementType, measurementUnit) %>%
  arrange(measurementType)

print(type_unit_pairs, n = Inf)

depth_lookup <- measurements %>%
  filter(grepl("depth|dyp|pressure", measurementType, ignore.case = TRUE)) %>%
  mutate(depth_val = readr::parse_number(as.character(measurementValue))) %>%
  filter(!is.na(depth_val)) %>%
  group_by(id) %>%
  summarise(depth_val = mean(depth_val, na.rm = TRUE), .groups = "drop")

# Find events with surface samples (<= 2m)
# We join by eventID to ensure we don't lose parameters that don't have their own depth row
gbif_historical_sf <- measurements %>%
  # 1. Use left_join so we don't lose parameters that lack a depth row
  left_join(depth_lookup, by = "id") %>% 
  
  # 2. Join to the event table for metadata (Lat/Lon/Date)
  left_join(events, by = c("id" = "eventID")) %>% 
  
  mutate(
    date = as.Date(eventDate),
    Year = year(date),
    # Use 0 as a fallback if depth is still NA (common for surface samples)
    depth = replace_na(depth_val, 0),
    
    parameter = case_when(
      grepl("nitrogen|Tot-N|Tot-N (µM)", measurementType, ignore.case = TRUE) ~ "Nitrogen",
      grepl("phosphorus|Tot-P|Tot-P(µM)", measurementType, ignore.case = TRUE) ~ "Phosphorus",
      grepl("secchi|siktedyp|Secchi depth", measurementType, ignore.case = TRUE) ~ "Secchi",
      grepl("temperature|temp", measurementType, ignore.case = TRUE) ~ "Temperature",
      grepl("salinity|salt|Salinity (PSU)", measurementType, ignore.case = TRUE) ~ "Salinity",
      grepl("NO3-N (µM)|nitrat", measurementType, ignore.case = TRUE) ~ "Nitrate",
      grepl("PO4-P(µM)|fosfat", measurementType, ignore.case = TRUE) ~ "Phosphate",
      grepl("Chl-a|Klorofyll", measurementType, ignore.case = TRUE) ~ "Chla",
      TRUE ~ "Other"
    ),
    value = as.numeric(measurementValue),
    source = "GBIF_Archive"
  ) %>%
  # 3. Clean up
  filter(parameter != "Other", !is.na(decimalLatitude)) %>%
  select(site = locality, date, Year, parameter, value, unit = measurementUnit, 
         source, depth, lat = decimalLatitude, lon = decimalLongitude) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


# --- 2. Load & Standardize Vannmiljø (Aggressive Parameter Matching) ---
param_mapping <- c(
  "Totalnitrogen" = "Nitrogen",
  "Totalfosfor"   = "Phosphorus",
  "Nitrat"         = "Nitrate",         # Typical Vannmiljø name
  "Nitrat + nitritt"= "Nitrate",         # Alternate name
  "Fosfat"         = "Phosphate",
  "Siktedyp"         = "Secchi",
  "Temperatur"     = "Temperature",
  "Salinitet"      = "Salinity",
  "Kjeldahl-nitrogen" = "orgN",
  "Løst organisk karbon (DOC)" = "DOC",
  "Totalt organisk karbon (TOC)" = "TOC",
  "Orto-fosfat" = "Phosphate",
  "Klorofyll a" = "Chla",
  "Fluoresensmåling av klorofyll a" = "Chla"
)

# Add your new file to this list
files <- c(
  "WaterRegistrationExport_f1.xlsx",
  "WaterRegistrationExport_f2.xlsx", # Your new file
  "WaterRegistrationExport_f3.xlsx"
)

# 2. Combine all files into one master dataframe
# map_df reads each file and stacks them on top of each other
df_combined <- files %>%
  set_names() %>% 
  map_df(~read_excel(.x), .id = "source_file")

vann_parameter_names <- df_combined %>%
  distinct(Parameter_navn) %>%
  arrange(Parameter_navn)

# Print the full list to the console
print(vann_parameter_names, n = Inf)

# --- Audit: Check Units paired with Parameters ---
# This is crucial for your standardization logic later
vann_unit_check <- df_combined %>%
  distinct(Parameter_navn, Enhet) %>%
  arrange(Parameter_navn)

print(vann_unit_check, n = Inf)

vann_sf <- df_combined %>%
  filter(Parameter_navn %in% names(param_mapping)) %>%
  # filter(Nedre_dyp <= 2 | is.na(Nedre_dyp)) %>%
  mutate(
    date = as.Date(as.POSIXct(Tid_provetak, format="%Y-%m-%d %H:%M:%S")),
    Year = year(date),
    parameter = param_mapping[Parameter_navn],
    source = "Vannmiljø",
    value = as.numeric(Verdi),
    depth = as.numeric(Nedre_dyp)
  ) %>%
  # Dynamic coordinate selection in case of naming variations
  rename(x = contains("Ost"), y = contains("Nord")) %>% 
  filter(!is.na(x)) %>%
  st_as_sf(coords = c("x", "y"), crs = 25833) %>%
  st_transform(4326)

# --- 3. Master Merge & Unit Standardization ---
oslo_final_sf <- bind_rows(gbif_historical_sf, vann_sf) %>%
  # Ensure everything is in 4326 for the Fjord_Part logic
  st_transform(4326) %>%
  mutate(
    # Fix coordinate extraction so it works for both sources
    lon_val = st_coordinates(.)[,1],
    lat_val = st_coordinates(.)[,2],
    
    # Standardize Units & Values
    needs_conv = grepl("µM|uM", unit, ignore.case = TRUE),
    value = case_when(
      parameter == "Phosphorus" & needs_conv ~ value * 30.97,
      parameter == "Nitrogen" & needs_conv ~ value * 14.01,
      TRUE ~ value
    ),
    unit = case_when(
      parameter == "Salinity" ~ "PSU",
      parameter %in% c("Nitrogen", "Phosphorus") & needs_conv ~ "µg/l",
      TRUE ~ unit
    )
  ) %>%
  # Season and Spatial Logic (Unified)
  mutate(
    Month = month(date),
    Season = case_when(
      Month %in% c(12, 1, 2) ~ "Vinter",
      Month %in% c(3, 4, 5)  ~ "Vår",
      Month %in% c(6, 7, 8)  ~ "Sommer",
      Month %in% c(9, 10, 11) ~ "Høst"
    ),
    Fjord_Part = case_when(
      lat_val > 59.50 & lon_val < 10.45 ~ "Drammensfjord",
      lat_val > 59.68 ~ "Inner Oslofjord",
      lon_val <= 10.60 ~ "Outer West",
      lon_val > 10.60 ~ "Outer East",
      TRUE ~ "Unclassified"
    )
  )

# --- CRITICAL CHECK: Did GBIF survive the merge? ---
table(oslo_final_sf$source)

# Temperal validation
ggplot(oslo_final_sf, aes(x = Year, fill = source)) +
  geom_histogram(binwidth = 1, position = "stack") +
  facet_wrap(~parameter, scales = "free_y") +
  theme_minimal() +
  labs(title = "Validation: Temporal Data Distribution", y = "Observations")

# Spatial validation
oslo_final_sf %>%
  st_drop_geometry() %>%
  count(Fjord_Part, parameter, source) %>%
  ggplot(aes(x = Fjord_Part, y = n, fill = source)) +
  geom_col() +
  facet_wrap(~parameter) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Validation: Spatial Distribution across Oslo Fjord")

# Data Loss Report
cat("--- FINAL DATA AUDIT ---\n")
cat("GBIF Archive points recovered: ", nrow(gbif_historical_sf), "\n")
cat("Vannmiljø points recovered:    ", nrow(vann_sf), "\n")
cat("Total Surface Samples:         ", nrow(oslo_final_sf), "\n")

# cluster sampling sites based on 1 km radius
vann_projected <- oslo_final_sf %>% 
  st_transform(25833)
coords <- st_coordinates(vann_projected)

# Calculate distance and cluster
# 'complete' linkage ensures the maximum distance within a cluster is < h
# 1. Project to metric for clustering
oslo_projected <- oslo_final_sf %>% st_transform(25833)
coords_for_db <- st_coordinates(oslo_projected)

# 2. Run DBSCAN (1000m radius)
# minPts = 1 ensures historical single-point stations aren't ignored
db_results <- dbscan(coords_for_db, eps = 1000, minPts = 1)
oslo_final_sf$site_cluster <- db_results$cluster

# 3. Create names (Handles cases where GBIF has 'site' and Vannmiljø has 'Vannlokalitetsnavn')
oslo_final_sf <- oslo_final_sf %>%
  group_by(site_cluster) %>%
  mutate(
    cluster_name = case_when(
      # 1. Try Vannmiljø names first
      any(!is.na(Vannlokalitetsnavn)) ~ {
        nm <- names(sort(table(Vannlokalitetsnavn, useNA = "no"), decreasing = TRUE))
        if(length(nm) > 0) nm[1] else NA_character_
      },
      # 2. Try GBIF site names second
      any(!is.na(site)) ~ {
        nm <- names(sort(table(site, useNA = "no"), decreasing = TRUE))
        if(length(nm) > 0) nm[1] else NA_character_
      },
      # 3. Fallback to Cluster ID
      TRUE ~ paste("Cluster", first(site_cluster))
    )
  ) %>%
  # Fill any remaining NAs with the Cluster ID just in case
  mutate(cluster_name = coalesce(cluster_name, paste("Cluster", site_cluster))) %>%
  ungroup()

# --- 2. Validation Map ---
# Calculate clean x-axis breaks
lon_range <- st_bbox(oslo_final_sf)
x_breaks <- seq(floor(lon_range["xmin"]*5)/5, ceiling(lon_range["xmax"]*5)/5, by = 0.2)

ggplot(oslo_final_sf) +
  # Detailed base map
  annotation_map_tile(type = "cartolight", zoom = 10, alpha = 0.9) + 
  
  # Map points
  layer_spatial(data = oslo_final_sf, aes(color = Fjord_Part), size = 1.5, alpha = 0.8) +
  
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = x_breaks) + 
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    # Margin at bottom for the scale and arrow
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 60, l = 10, unit = "pt")
  ) +
  
  labs(
    title = "Validation: Geographic Fjord Regions",
    subtitle = "Confirmed regions using standardized master data",
    x = "Longitude (°E)", 
    y = "Latitude (°N)"
  ) +
  
  # Map elements outside the frame
  annotation_scale(
    location = "br", 
    pad_y = unit(-0.45, "in"),
    text_cex = 0.7
  ) +
  annotation_north_arrow(
    location = "br", 
    pad_y = unit(-0.75, "in"),
    style = north_arrow_orienteering(text_size = 8),
    height = unit(0.8, "cm"),
    width = unit(0.8, "cm")
  ) +
  coord_sf(clip = "off")

# Subset to data from 2000 and onwards
oslo_post_2000 <- oslo_final_sf %>%
  filter(Year >= 2000)

# Optional: Verify the date range
summary(oslo_post_2000$Year)

oslo_final_sf <- oslo_final_sf %>%
  mutate(
    # Create the weekly bin (Monday-start)
    Week_Date = floor_date(as.Date(date), unit = "week", week_start = 1),
    # Ensure Year is numeric
    Year = year(date)
  )

df_ecosystem_3d <- oslo_final_sf %>%
  st_drop_geometry() %>%
  mutate(
    # 7-day temporal window
    Week_Date = floor_date(date, "week"),
    # Vertical Bins for the manuscript
    Depth_Bin = case_when(
      depth <= 5 ~ "0-5m",
      depth > 5 & depth <= 20 ~ "5-20m (DCM Zone)",
      depth > 20 & depth <= 50 ~ "20-50m",
      depth > 50 ~ ">50m",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Aggregate to reduce memory before pivoting
  group_by(Year, Week_Date, Fjord_Part, site_cluster, Depth_Bin, parameter) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value)

# --- 5. Tracking the DCM ---
# Find the depth of maximum chlorophyll for every week/cluster
df_dcm_track <- df_ecosystem_3d %>%
  filter(!is.na(Chla)) %>%
  group_by(Year, Week_Date, site_cluster) %>%
  slice_max(Chla, n = 1, with_ties = FALSE) %>%
  select(Year, Week_Date, site_cluster, DCM_Depth = Depth_Bin, DCM_Value = Chla)
