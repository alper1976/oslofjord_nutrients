library(tidyverse)

# Check 1: Are the coordinates in degrees or millions?
summary(oslo_final_sf$lat_val) 

# Check 2: How many rows exist per region before filtering?
table(oslo_final_sf$Fjord_Part, useNA = "always")

# Check 3: Is it a parameter issue?
table(oslo_final_sf$Fjord_Part, oslo_final_sf$parameter)

# How much data do you actually have?
oslo_final_sf %>%
  st_drop_geometry() %>%
  group_by(Fjord_Part, parameter) %>%
  summarise(unique_years = n_distinct(Year), total_samples = n()) %>%
  arrange(unique_years)

# Evaluation: Where is the data?
coverage_check <- oslo_final_sf %>%
  st_drop_geometry() %>%
  mutate(Decade = (Year %/% 10) * 10) %>%
  count(Fjord_Part, Decade, parameter)

ggplot(coverage_check, aes(x = Decade, y = Fjord_Part, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "mako", trans = "log10", name = "Log(Samples)") +
  facet_wrap(~parameter) +
  theme_minimal() +
  labs(
    title = "Data Density Evaluation: Fjord_Part vs. Time",
    subtitle = "White space indicates total lack of data in that region/decade",
    x = "Decade",
    y = "Fjord Region"
  )

#Evaluate depth distribution
table(oslo_final_sf$source, oslo_final_sf$depth > 20)
# 1. Plot the full depth distribution
oslo_final_sf %>%
  filter(depth <= 400) %>% # Filter out any extreme deep outliers
  ggplot(aes(x = depth, fill = parameter)) +
  
  # Use a 5-meter bin for the water column
  geom_histogram(binwidth = 5, alpha = 0.8, show.legend = FALSE) +
  
  # Log scale for the count so we can see the deep measurements
  scale_y_log10() +
  
  # Reverse X-axis (0 at top, 400 at bottom) and enforce the limits
  scale_x_reverse(limits = c(400, 0)) + 
  coord_flip() +
  
  facet_wrap(~parameter, scales = "free_x", ncol = 4) +
  theme_minimal() +
  labs(
    title = "Water Column Distribution: 0m to 400m",
    subtitle = "Y-axis (Observations) is Log10 scaled to highlight deep data points",
    x = "Depth (meters)",
    y = "Log10(Count of Observations)"
  ) +
  theme(strip.text = element_text(face = "bold"))
# 2. Print a quick numeric check of your maximum depths per parameter
oslo_final_sf %>%
  st_drop_geometry() %>%
  group_by(parameter) %>%
  summarise(
    Max_Depth_m = max(depth, na.rm = TRUE),
    Observations_below_100m = sum(depth > 100, na.rm = TRUE),
    .groups = "drop"
  )

## Get top 5 sites when it comes to data points.
top_5_sites <- oslo_final_sf %>%
  st_drop_geometry() %>%
  count(Fjord_Part, cluster_name, sort = TRUE) %>%
  group_by(Fjord_Part) %>%
  slice_max(n, n = 5)

print(top_5_sites)