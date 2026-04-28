# ==============================================================================
# SATELLITE INTEGRATION MASTER SCRIPT
# Objective: Mine NASA MODIS data, extract bimodal peaks via GAMs, validate 
# historical in situ phenology shifts, and compare multi-decadal trends.
# ==============================================================================

library(rerddap)
library(tidyverse)
library(sf)
library(lubridate)
library(mgcv)
library(broom)

# ---------------------------------------------------------
# 1. ERDDAP Query: Download 8-Day MODIS-Aqua Composites
# ---------------------------------------------------------
cat("\n[1/7] Contacting NASA ERDDAP for 8-Day MODIS Chla...\n")

# Let's mine the last 20 years to match your modern data
start_date <- "2003-02-16" 
end_date   <- "2022-04-16"

sat_8day <- griddap(
  datasetx = "erdMH1chla8day",
  time = c(start_date, end_date),
  latitude = c(59.0, 59.9),
  longitude = c(10.2, 10.8),
  fields = "chlorophyll"
)

# ---------------------------------------------------------
# 2. Spatial Matching: Link Pixels to Fjord Regions
# ---------------------------------------------------------
cat("\n[2/7] Matching satellite pixels to established Fjord Regions...\n")

# Extract unique coordinates and regions from your boat data
cluster_coords <- oslo_final_sf %>%
  bind_cols(as_tibble(st_coordinates(.))) %>% 
  st_drop_geometry() %>%
  rename(cluster_lon = X, cluster_lat = Y) %>% 
  distinct(site_cluster, cluster_lon, cluster_lat, Fjord_Part) 
# Inner Oslofjord is perfectly preserved here!

# Clean the raw ERDDAP data and cross-reference with our regions
df_sat_matched <- sat_8day$data %>%
  drop_na(chlorophyll) %>%
  rename(sat_lon = longitude, sat_lat = latitude) %>%
  mutate(
    # Stripping hidden matrices that crash the console
    sat_lon = as.numeric(sat_lon),
    sat_lat = as.numeric(sat_lat),
    chlorophyll = as.numeric(chlorophyll),
    date = as.Date(time),
    Year = year(date),
    DOY = yday(date),
    Month = month(date)
  ) %>%
  # Join and find the closest pixel to our clusters
  cross_join(cluster_coords) %>%
  mutate(dist = sqrt((cluster_lon - sat_lon)^2 + (cluster_lat - sat_lat)^2)) %>%
  group_by(site_cluster, date) %>%
  slice_min(dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # Calculate the regional average for each 8-day window
  group_by(Fjord_Part, Year, DOY, Month) %>%
  summarise(Mean_Chla = mean(chlorophyll, na.rm = TRUE), .groups = "drop")


# ---------------------------------------------------------
# 3. Fit GAMs to Extract Bimodal Peaks (Spring vs Summer)
# ---------------------------------------------------------
cat("\n[3/7] Fitting regional GAMs to extract phenological peaks...\n")

peak_comparison <- tibble(
  Fjord_Part = character(),
  Year = integer(), 
  Spring_Peak_DOY = numeric(), 
  Spring_Peak_Chla = numeric(),
  Summer_Peak_DOY = numeric(),
  Summer_Peak_Chla = numeric()
)

# Loop through each Region and Year
for(region in unique(df_sat_matched$Fjord_Part)) {
  for(yr in unique(df_sat_matched$Year)) {
    
    yr_data <- df_sat_matched %>% filter(Fjord_Part == region, Year == yr)
    
    # Require at least 15 clear 8-day composite readings across the year
    if(nrow(yr_data) >= 15) { 
      
      fit <- gam(Mean_Chla ~ s(DOY, bs = "cs", k = 10), data = yr_data)
      daily_seq <- data.frame(DOY = seq(1, 365, by = 1))
      
      # Predict values for the entire year
      predicted_chla <- as.numeric(predict(fit, newdata = daily_seq)) 
      daily_predictions <- tibble(DOY = daily_seq$DOY, Pred_Chla = predicted_chla)
      
      # Define biological windows
      spring_window <- daily_predictions %>% filter(DOY >= 32 & DOY <= 151)   
      summer_window <- daily_predictions %>% filter(DOY >= 166 & DOY <= 304)  
      
      # Extract maximums 
      spring_max <- spring_window %>% slice_max(Pred_Chla, n = 1, with_ties = FALSE)
      summer_max <- summer_window %>% slice_max(Pred_Chla, n = 1, with_ties = FALSE)
      
      peak_comparison <- bind_rows(peak_comparison, tibble(
        Fjord_Part = region,
        Year = yr,
        Spring_Peak_DOY = spring_max$DOY,
        Spring_Peak_Chla = spring_max$Pred_Chla,
        Summer_Peak_DOY = summer_max$DOY,
        Summer_Peak_Chla = summer_max$Pred_Chla
      ))
    }
  }
}

peak_comparison <- peak_comparison %>% mutate(Peak_Ratio = Spring_Peak_Chla / Summer_Peak_Chla)


# ---------------------------------------------------------
# 4. Validation Plot: Boat vs Satellite Phenology
# ---------------------------------------------------------
cat("\n[4/7] Generating Validation plots...\n")

df_validation <- df_bloom_timing %>%
  select(Fjord_Part, Year, InSitu_Bloom_DOY = Peak_DOY) %>%
  inner_join(peak_comparison %>% select(Fjord_Part, Year, Sat_Bloom_DOY = Spring_Peak_DOY), 
             by = c("Fjord_Part", "Year"))

p_sat_validation <- ggplot(df_validation) +
  geom_point(aes(x = Year, y = InSitu_Bloom_DOY), color = "#377eb8", alpha = 0.5, size = 2) +
  geom_smooth(aes(x = Year, y = InSitu_Bloom_DOY, color = "In Situ (Boat)"), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 3), se = FALSE, linewidth = 1.2) +
  
  geom_point(aes(x = Year, y = Sat_Bloom_DOY), color = "#e41a1c", alpha = 0.5, size = 2, shape = 17) +
  geom_smooth(aes(x = Year, y = Sat_Bloom_DOY, color = "Satellite (Ocean Color)"), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 3), se = FALSE, linewidth = 1.2) +
  
  geom_hline(yintercept = c(32, 60, 91, 121), linetype = "dotted", color = "gray60") +
  facet_wrap(~Fjord_Part) +
  theme_bw() +
  scale_color_manual(values = c("In Situ (Boat)" = "#377eb8", "Satellite (Ocean Color)" = "#e41a1c")) +
  labs(
    title = "Validation of Spring Bloom Phenology",
    subtitle = "Comparing historical in situ data with modern satellite extraction",
    y = "Day of Year (Julian Day)", x = "Year", color = "Observation Method"
  ) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

print(p_sat_validation)


# ---------------------------------------------------------
# 5. Plot: Satellite Bimodal Regime Shift
# ---------------------------------------------------------
cat("\n[5/7] Evaluating Regime Shifts...\n")

p_sat_regime <- ggplot(peak_comparison, aes(x = Year, y = Peak_Ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_point(aes(color = Peak_Ratio > 1), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), color = "black") +
  
  facet_wrap(~Fjord_Part) + theme_bw() +
  scale_color_manual(values = c("TRUE" = "#41ab5d", "FALSE" = "#ef3b2c")) +
  labs(
    title = "Satellite-Derived Regime Shift (2003-2023)",
    subtitle = "Values > 1 indicate Spring dominance. Values < 1 indicate Summer dominance.",
    y = "Ratio (Spring Peak / Summer Peak)", x = "Year"
  ) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

print(p_sat_regime)


# ---------------------------------------------------------
# 6. Satellite vs. In Situ Trend Comparison
# ---------------------------------------------------------
cat("\n[6/7] Generating Satellite vs. In Situ Trend Comparison...\n")

# A. Prepare the In Situ Surface Data (Annual Means)
df_insitu_surface <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(grepl("Chl", parameter, ignore.case = TRUE), depth <= 5) %>%
  mutate(
    Year = year(date), 
    value = as.numeric(value),
    Fjord_Part = as.character(Fjord_Part) # Force character to ensure safe facetting
  ) %>%
  filter(is.finite(value), value > 0, Year >= 2003) %>% 
  group_by(Fjord_Part, Year) %>%
  summarise(InSitu_Chla = mean(value, na.rm = TRUE), .groups = "drop")

# B. Prepare the Satellite Data
df_sat <- df_sat_matched %>%
  mutate(Fjord_Part = as.character(Fjord_Part)) %>% # Force character matching
  group_by(Fjord_Part, Year) %>%
  summarise(Sat_Chla = mean(Mean_Chla, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(Sat_Chla))

# DIAGNOSTIC CHECK: Prove the satellite data exists in the console
cat("\n--- CHECK: First 5 rows of df_sat (Satellite Data) ---\n")
print(head(df_sat))

# C. Define the new scaling factor
# Scaled to 5x as requested
scale_factor <- 5 

# Calculate a dynamic upper limit so data is NEVER accidentally clipped
max_sat <- max(df_sat$Sat_Chla, na.rm = TRUE)
max_insitu <- max(df_insitu_surface$InSitu_Chla * scale_factor, na.rm = TRUE)
dynamic_upper_limit <- max(10, max_sat, max_insitu) * 1.1 # Expands automatically based on the 5x factor

# D. Plot Independent Dataframes on Dual Axes
p_sat_vs_insitu <- ggplot() +
  
  # Layer 1: In Situ Data (Scaled 5x)
  geom_point(data = df_insitu_surface, aes(x = Year, y = InSitu_Chla * scale_factor), 
             color = "#377eb8", alpha = 0.4, size = 1.5) +
  geom_smooth(data = df_insitu_surface, aes(x = Year, y = InSitu_Chla * scale_factor), 
              method = "lm", color = "#377eb8", linetype = "dashed", linewidth = 1.2, se = FALSE) +
  
  # Layer 2: Satellite Data (Unscaled)
  geom_point(data = df_sat, aes(x = Year, y = Sat_Chla), 
             color = "#e41a1c", alpha = 0.6, size = 2) +
  geom_smooth(data = df_sat, aes(x = Year, y = Sat_Chla), 
              method = "lm", color = "#e41a1c", linewidth = 1.2, se = TRUE, alpha = 0.2) +
  
  # Apply the dynamic limit so nothing is deleted
  coord_cartesian(ylim = c(0, dynamic_upper_limit)) +
  
  # Dual Y-Axis Setup (Calculates the Right Axis correctly based on the scale_factor)
  scale_y_continuous(
    name = "Satellite Chla (mg/m³)", 
    sec.axis = sec_axis(~ . / scale_factor, name = "In Situ Chla (µg/L)")
  ) +
  
  facet_wrap(~ Fjord_Part) + 
  
  theme_bw() +
  labs(
    title = "Relative Surface Chlorophyll Trends (2003-2023)",
    subtitle = paste0("Red Solid = Satellite (MODIS) | Blue Dashed = In Situ (Boat, Scaled ", scale_factor, "x)"),
    x = "Year"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.title.y.left = element_text(color = "#e41a1c", face = "bold", margin = ggplot2::margin(r = 10)),
    axis.title.y.right = element_text(color = "#377eb8", face = "bold", margin = ggplot2::margin(l = 10)),
    panel.grid.minor = element_blank()
  )

print(p_sat_vs_insitu)

# ---------------------------------------------------------
# 7. Print/Export Files to System
# ---------------------------------------------------------
cat("\n[7/7] Printing plots and exporting data files...\n")

# Create directories if they don't exist
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Export Plots as PDFs
ggsave("outputs/figures/Fig_Sat_Validation.pdf", plot = p_sat_validation, width = 10, height = 8)
ggsave("outputs/figures/Fig_Sat_Regime_Shift.pdf", plot = p_sat_regime, width = 10, height = 8)
ggsave("outputs/figures/Supplementary_Figure_8_Satellite_A4.pdf", plot = p_sat_vs_insitu, width = 12, height = 8)

# Export processed datasets
saveRDS(df_sat_matched, "data/processed/satellite_cleaned.rds")
write_csv(peak_comparison, "data/processed/satellite_phenology_peaks.csv")

cat("\nDone! All files successfully printed to 'outputs/figures/' and 'data/processed/'.\n")

