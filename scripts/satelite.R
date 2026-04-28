# ==============================================================================
# SATELLITE INTEGRATION MASTER SCRIPT
# Objective: Mine NASA MODIS data, extract bimodal peaks via GAMs, and 
# validate historical in situ phenology shifts.
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
cat("\n[1/4] Contacting NASA ERDDAP for 8-Day MODIS Chla...\n")


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
cat("\n[2/4] Matching satellite pixels to established Fjord Regions...\n")

# Extract unique coordinates and regions from your boat data
cluster_coords <- oslo_final_sf %>%
  bind_cols(as_tibble(st_coordinates(.))) %>% 
  st_drop_geometry() %>%
  rename(cluster_lon = X, cluster_lat = Y) %>% 
  distinct(site_cluster, cluster_lon, cluster_lat, Fjord_Part) %>%
  filter(Fjord_Part != "Inner Oslofjord") # Match the boat data constraints

# Clean the raw ERDDAP data and cross-reference with our regions
df_sat_matched <- sat_8day$data %>%
  drop_na(chlorophyll) %>%
  rename(sat_lon = longitude, sat_lat = latitude) %>%
  mutate(
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
cat("\n[3/4] Fitting regional GAMs to extract phenological peaks...\n")

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
    
    # Require at least 15 clear 8-day composite readings across the year to fit a double-hump curve
    if(nrow(yr_data) >= 15) { 
      
      # Use k=10 so the curve is flexible enough to bend up and down twice
      fit <- gam(Mean_Chla ~ s(DOY, bs = "cs", k = 10), data = yr_data)
      
      # Predict values for the entire year
      daily_seq <- data.frame(DOY = seq(1, 365, by = 1))
      predicted_chla <- predict(fit, newdata = daily_seq)
      daily_predictions <- tibble(DOY = daily_seq$DOY, Pred_Chla = predicted_chla)
      
      # Define biological windows
      spring_window <- daily_predictions %>% filter(DOY >= 32 & DOY <= 151)   # Feb 1 - May 31
      summer_window <- daily_predictions %>% filter(DOY >= 166 & DOY <= 304)  # Jun 15 - Oct 31
      
      # Extract maximums
      spring_max <- spring_window %>% slice_max(Pred_Chla, n = 1)
      summer_max <- summer_window %>% slice_max(Pred_Chla, n = 1)
      
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

# Calculate Regime Ratio (> 1 means Spring dominates)
peak_comparison <- peak_comparison %>%
  mutate(Peak_Ratio = Spring_Peak_Chla / Summer_Peak_Chla)


# ---------------------------------------------------------
# 4. Validation Plot: Boat vs Satellite Phenology
# ---------------------------------------------------------
cat("\n[4/4] Generating Validation plots...\n")

# NOTE: This assumes 'df_bloom_timing' is already in your environment from the In Situ script!
df_validation <- df_bloom_timing %>%
  select(Fjord_Part, Year, InSitu_Bloom_DOY = Bloom_DOY) %>%
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
    y = "Day of Year (Julian Day)",
    x = "Year",
    color = "Observation Method"
  ) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

print(p_sat_validation)

# Calculate Validation Correlation
validation_stats <- df_validation %>%
  group_by(Fjord_Part) %>%
  summarise(
    Pearson_r = cor(InSitu_Bloom_DOY, Sat_Bloom_DOY, method = "pearson", use = "complete.obs"),
    .groups = "drop"
  )
cat("\n--- PEARSON CORRELATION: SATELLITE VS BOAT ---\n")
print(validation_stats)


# ---------------------------------------------------------
# 5. Plot: Satellite Bimodal Regime Shift
# ---------------------------------------------------------
p_sat_regime <- ggplot(peak_comparison, aes(x = Year, y = Peak_Ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_point(aes(color = Peak_Ratio > 1), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), color = "black") +
  
  facet_wrap(~Fjord_Part) +
  theme_bw() +
  scale_color_manual(values = c("TRUE" = "#41ab5d", "FALSE" = "#ef3b2c")) +
  labs(
    title = "Satellite-Derived Regime Shift (2003-2023)",
    subtitle = "Values > 1 indicate Spring dominance. Values < 1 indicate Summer dominance.",
    y = "Ratio (Spring Peak / Summer Peak)",
    x = "Year"
  ) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

print(p_sat_regime)

# Extract significance for the Satellite regime shift
cat("\n--- SATELLITE REGIME SHIFT TRENDS ---\n")
sat_regime_stats <- peak_comparison %>%
  group_by(Fjord_Part) %>%
  reframe(tidy(lm(Peak_Ratio ~ Year, data = cur_data()))) %>%
  filter(term == "Year") %>%
  mutate(Significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  select(Fjord_Part, Slope = estimate, p_value = p.value, Significant)

print(sat_regime_stats)

# ---------------------------------------------------------
# 6. Biomass Trend Validation: Dual-Scale Annual Means
# ---------------------------------------------------------
cat("\n[6/4] Evaluating multi-method Chla decline with dual-scale normalization...\n")

# 1. Standardize and Combine Data
# Using the dfs generated in the previous step (df_insitu_annual and df_sat_annual)
df_trend_comp <- bind_rows(df_insitu_annual, df_sat_annual)

# 2. Plot with Dual Y-Axis Scaling
# We multiply In Situ by 10 to map it onto the Satellite scale for visual comparison
p_chla_dual_scale <- ggplot(df_trend_comp, aes(x = Year, y = Mean_Chla, color = Method)) +
  # In Situ points and line (Scaled up 10x for the plot)
  geom_point(data = filter(df_trend_comp, Method == "In Situ (Boat)"), 
             aes(y = Mean_Chla * 10), alpha = 0.5) +
  geom_smooth(data = filter(df_trend_comp, Method == "In Situ (Boat)"), 
              aes(y = Mean_Chla * 10), method = "lm", se = TRUE, linetype = "dashed") +
  
  # Satellite points and line (Original Scale)
  geom_point(data = filter(df_trend_comp, Method == "Satellite (MODIS)"), 
             alpha = 0.5) +
  geom_smooth(data = filter(df_trend_comp, Method == "Satellite (MODIS)"), 
              method = "lm", se = TRUE) +
  
  facet_wrap(~Fjord_Part, scales = "free_y") +
  
  # Create the Dual Axis
  scale_y_continuous(
    name = "Satellite Chla (mg/m³)",
    sec.axis = sec_axis(~ . / 10, name = "In Situ Chla (µg/L)")
  ) +
  
  scale_color_manual(values = c("In Situ (Boat)" = "#377eb8", "Satellite (MODIS)" = "#e41a1c")) +
  theme_bw() +
  labs(
    title = "Relative Surface Chlorophyll Trends (2003-2023)",
    subtitle = "In Situ values scaled 10x to facilitate trend comparison",
    x = "Year"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(), # Removes "Method" title
    strip.text = element_text(face = "bold"),
    axis.title.y.right = element_text(color = "#377eb8"),
    axis.title.y.left = element_text(color = "#e41a1c")
  )

print(p_chla_dual_scale)

# 3. Statistical Comparison of Slopes
# We compare the percentage change per decade to normalize the units
stats_comp <- df_trend_comp %>%
  group_by(Fjord_Part, Method) %>%
  do(tidy(lm(Mean_Chla ~ Year, data = .))) %>%
  filter(term == "Year") %>%
  group_by(Fjord_Part) %>%
  mutate(
    Relative_Decadal_Change = (estimate * 10) / mean(df_trend_comp$Mean_Chla[df_trend_comp$Fjord_Part == first(Fjord_Part) & df_trend_comp$Method == first(Method)]) * 100
  )

cat("\n--- RELATIVE BIOMASS TRENDS (% Change per Decade) ---\n")
print(stats_comp %>% select(Fjord_Part, Method, p_value = p.value, Rel_Change_Pct = Relative_Decadal_Change))
