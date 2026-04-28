
library(ggplot2)
library(broom)
library(purrr)

# ==============================================================================
# PYCNOCLINE: Prepare data
# ==============================================================================
# Salinity Pycnocline Script (Depth & Strength Filtered)
df_pycnocline <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(parameter == "Salinity", !is.na(value), !is.na(depth)) %>%
  group_by(Fjord_Part, Year, Season, Week_Date, site_cluster, depth) %>%
  summarise(Salinity = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  group_by(Fjord_Part, Year, Season, Week_Date, site_cluster) %>%
  arrange(depth, .by_group = TRUE) %>%
  filter(n_distinct(depth) >= 2) %>%
  mutate(
    delta_S = Salinity - lag(Salinity),
    delta_Z = depth - lag(depth),
    gradient = delta_S / delta_Z
  ) %>%
  
  # Physical Filter: Require at least 0.5m between readings
  filter(!is.na(gradient), is.finite(gradient), delta_Z >= 0.5) %>%
  
  slice_max(gradient, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # Statistical Filter (IQR) for BOTH Strength and Depth
  group_by(Fjord_Part) %>%
  mutate(
    # Calculate bounds for Strength (Gradient)
    Q1_grad = quantile(gradient, 0.25, na.rm = TRUE),
    Q3_grad = quantile(gradient, 0.75, na.rm = TRUE),
    IQR_grad = Q3_grad - Q1_grad,
    Upper_grad = Q3_grad + (3 * IQR_grad),
    
    # Calculate bounds for Depth
    Q1_depth = quantile(depth, 0.25, na.rm = TRUE),
    Q3_depth = quantile(depth, 0.75, na.rm = TRUE),
    IQR_depth = Q3_depth - Q1_depth,
    Lower_depth = Q1_depth - (3 * IQR_depth), 
    Upper_depth = Q3_depth + (3 * IQR_depth)
  ) %>%
  # Apply both filters (only keep natural gradients at natural depths)
  filter(
    gradient <= Upper_grad, 
    depth >= Lower_depth, 
    depth <= Upper_depth
  ) %>%
  
  select(Fjord_Part, Year, Season, Week_Date, site_cluster, 
         Pycnocline_Depth = depth, Max_Gradient = gradient) %>%
  ungroup()

# True Density Pycnocline Script (Depth & Strength Filtered)

df_true_pycnocline <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(parameter %in% c("Salinity", "Temperature"), !is.na(value), !is.na(depth)) %>%
  group_by(Fjord_Part, Year, Season, Week_Date, site_cluster, depth, parameter) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  filter(!is.na(Salinity), !is.na(Temperature)) %>%
  mutate(
    Density = swSigmaT(salinity = Salinity, temperature = Temperature, pressure = depth)
  ) %>%
  group_by(Fjord_Part, Year, Season, Week_Date, site_cluster) %>%
  arrange(depth, .by_group = TRUE) %>%
  filter(n_distinct(depth) >= 2) %>%
  mutate(
    delta_Rho = Density - lag(Density),
    delta_Z = depth - lag(depth),
    Density_Gradient = delta_Rho / delta_Z
  ) %>%
  
  # Physical Filter
  filter(!is.na(Density_Gradient), is.finite(Density_Gradient), delta_Z >= 0.5) %>%
  
  slice_max(Density_Gradient, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # Statistical Filter (IQR) for BOTH Strength and Depth
  group_by(Fjord_Part) %>%
  mutate(
    # Bounds for Density Strength
    Q1_grad = quantile(Density_Gradient, 0.25, na.rm = TRUE),
    Q3_grad = quantile(Density_Gradient, 0.75, na.rm = TRUE),
    IQR_grad = Q3_grad - Q1_grad,
    Upper_grad = Q3_grad + (3 * IQR_grad),
    
    # Bounds for Pycnocline Depth
    Q1_depth = quantile(depth, 0.25, na.rm = TRUE),
    Q3_depth = quantile(depth, 0.75, na.rm = TRUE),
    IQR_depth = Q3_depth - Q1_depth,
    Lower_depth = Q1_depth - (3 * IQR_depth),
    Upper_depth = Q3_depth + (3 * IQR_depth)
  ) %>%
  # Apply filters
  filter(
    Density_Gradient <= Upper_grad,
    depth >= Lower_depth,
    depth <= Upper_depth
  ) %>%
  
  select(Fjord_Part, Year, Season, Week_Date, site_cluster, 
         Pycnocline_Depth = depth, Pycnocline_Strength = Density_Gradient) %>%
  ungroup()

# ==============================================================================
# PYCNOCLINE DEPTH: TRENDS AND BOXPLOTS
# ==============================================================================

# --- A. Prepare Data for Trends (Require 5 years of data per site) ---
df_pyc_trends <- df_pycnocline %>%
  group_by(site_cluster) %>%
  filter(n_distinct(Year) >= 5) %>%
  ungroup()

# --- B. Plot 1: Site-Specific Pycnocline Trajectories ---
p_pyc_depth_trend <- ggplot(df_pyc_trends, aes(x = Year, y = Pycnocline_Depth, group = site_cluster)) +
  # Faded regression line for EACH site
  geom_smooth(method = "lm", se = FALSE, alpha = 0.4, aes(color = Fjord_Part), linewidth = 0.5) +
  # Overall regional average
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 1.5) +
  # Reverse Y-axis so surface (0) is at the top
  scale_y_reverse() + 
  facet_wrap(~Fjord_Part) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Site-Specific Pycnocline Depth Trajectories",
    subtitle = "Filtered for artifacts. Faded lines = Individual 1km sites. Black line = Regional average.",
    y = "Pycnocline Depth (meters)",
    x = "Year"
  ) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

print(p_pyc_depth_trend)

# --- C. Calculate Linear Models (Slopes) per Site ---
pyc_slopes <- df_pyc_trends %>%
  group_by(Fjord_Part, site_cluster) %>%
  nest() %>% 
  mutate(
    model = map(data, ~lm(Pycnocline_Depth ~ Year, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Year") %>%
  mutate(
    # Change per decade
    Decadal_Change = estimate * 10,
    Significant = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not Significant")
  )

# --- D. Plot 2: Boxplot Summarizing Decadal Rates ---
p_pyc_boxplot <- ggplot(pyc_slopes, aes(x = Fjord_Part, y = Decadal_Change)) +
  # Zero line (No change in depth)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  # Boxplot for regional spread
  geom_boxplot(fill = "lightblue", color = "gray40", outlier.shape = NA, alpha = 0.5) +
  # Jittered points for individual sites
  geom_jitter(aes(color = Significant), width = 0.2, size = 2.5, alpha = 0.8) +
  
  # Reverse Y-axis so "Up" means "Getting Shallower" (negative numbers)
  scale_y_reverse() +
  
  theme_minimal() +
  scale_color_manual(values = c("Not Significant" = "gray70", "Significant (p < 0.05)" = "#005a32")) +
  labs(
    title = "Rate of Pycnocline Shift by Site Cluster",
    subtitle = "Points ABOVE the dashed line mean the pycnocline is getting SHALLOWER.",
    y = "Change in Depth (Meters per Decade)",
    x = "Fjord Region",
    color = "Statistical Significance"
  ) +
  theme(axis.text.x = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(p_pyc_boxplot)

# ==============================================================================
# PYCNOCLINE STRENGTH: TRENDS AND BOXPLOTS
# ==============================================================================

# --- A. Prepare Data for Trends (Require 5 years of data per site) ---
df_strength_trends <- df_true_pycnocline %>%
  group_by(site_cluster) %>%
  filter(n_distinct(Year) >= 5) %>%
  ungroup()

# --- B. Plot 1: Site-Specific Pycnocline Strength Trajectories ---
p_strength_trend <- ggplot(df_strength_trends, aes(x = Year, y = Pycnocline_Strength, group = site_cluster)) +
  # Faded regression line for EACH site
  geom_smooth(method = "lm", se = FALSE, alpha = 0.4, aes(color = Fjord_Part), linewidth = 0.5) +
  # Overall regional average
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 1.5) +
  
  # Allow Y-axes to vary, as inner fjords often have much stronger density gradients than open coasts
  facet_wrap(~Fjord_Part, scales = "free_y") +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Site-Specific Pycnocline Strength Trajectories",
    subtitle = "Faded lines = Individual 1km sites. Black line = Regional average.",
    y = "Density Gradient (kg/m³ per meter)",
    x = "Year"
  ) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

print(p_strength_trend)


# --- C. Calculate Linear Models (Slopes) per Site ---
strength_slopes <- df_strength_trends %>%
  group_by(Fjord_Part, site_cluster) %>%
  nest() %>% 
  mutate(
    model = map(data, ~lm(Pycnocline_Strength ~ Year, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Year") %>%
  mutate(
    # Change per decade
    Decadal_Change = estimate * 10,
    Significant = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not Significant")
  )


# --- D. Plot 2: Boxplot Summarizing Decadal Rates ---
p_strength_boxplot <- ggplot(strength_slopes, aes(x = Fjord_Part, y = Decadal_Change)) +
  # Zero line (No change in stratification strength)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # Boxplot for regional spread (Using a warmer color to distinguish from the depth plot)
  geom_boxplot(fill = "lightcoral", color = "gray40", outlier.shape = NA, alpha = 0.5) +
  
  # Jittered points for individual sites
  geom_jitter(aes(color = Significant), width = 0.2, size = 2.5, alpha = 0.8) +
  
  theme_minimal() +
  # Custom colors: Gray for non-significant, deep red for significant
  scale_color_manual(values = c("Not Significant" = "gray70", "Significant (p < 0.05)" = "#b30000")) +
  labs(
    title = "Rate of Pycnocline Strength Shift by Site Cluster",
    subtitle = "Points ABOVE the dashed line mean the density barrier is getting STRONGER.",
    y = "Change in Strength (kg/m³ per m, per Decade)",
    x = "Fjord Region",
    color = "Statistical Significance"
  ) +
  theme(axis.text.x = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(p_strength_boxplot)






