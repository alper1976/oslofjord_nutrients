library(tidyverse)
library(ggplot2)
library(lubridate) # Needed for floor_date()
library(broom)     # Needed for tidy()
library(sf)        # Needed for st_drop_geometry()

# ==============================================================================
# MODULE 2: DATA AVAILABILITY HEATMAP (Updated for all 3 Ratios)
# Objective: Visualize temporal and spatial gaps in Total vs. Inorganic data.
# ==============================================================================

# 1. Create wide data for paired comparisons
df_ratios <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(parameter %in% c("Nitrogen", "Phosphorus", "Nitrate", "Phosphate")) %>%
  mutate(Year_Month = floor_date(date, "month")) %>%
  group_by(Fjord_Part, site_cluster, Year_Month, parameter) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  rename(Cluster_ID = site_cluster)

# 2. Identify existing pairs for ALL THREE ratios
availability_comparison <- df_ratios %>%
  mutate(
    Month = factor(month(Year_Month), levels = 1:12, labels = month.abb),
    Year = year(Year_Month),
    
    # 1. Total N and Total P exist
    Has_TN_Ratio      = if(all(c("Nitrogen", "Phosphorus") %in% names(.))) !is.na(Nitrogen) & !is.na(Phosphorus) else FALSE,
    # 2. Nitrate and Total P exist
    Has_NO3_Ratio     = if(all(c("Nitrate", "Phosphorus") %in% names(.))) !is.na(Nitrate)  & !is.na(Phosphorus) else FALSE,
    # 3. Nitrate and Phosphate exist (Purely Inorganic)
    Has_NO3_PO4_Ratio = if(all(c("Nitrate", "Phosphate") %in% names(.))) !is.na(Nitrate) & !is.na(Phosphate) else FALSE
  ) %>%
  # Pivot all three checks into a long format
  pivot_longer(
    cols = c(Has_TN_Ratio, Has_NO3_Ratio, Has_NO3_PO4_Ratio), 
    names_to = "Ratio_Type", 
    values_to = "Exists"
  ) %>%
  filter(Exists == TRUE) %>%
  group_by(Fjord_Part, Year, Month, Ratio_Type) %>%
  summarise(Cluster_Count = n_distinct(Cluster_ID), .groups = "drop") %>%
  
  # Clean up names and set the factor order to match the other plots perfectly
  mutate(
    Ratio_Type = recode(Ratio_Type,
                        "Has_TN_Ratio" = "Total N : Total P",
                        "Has_NO3_Ratio" = "Nitrate : Total P",
                        "Has_NO3_PO4_Ratio" = "Inorganic NO3 : PO4"),
    Ratio_Type = factor(Ratio_Type, levels = c("Total N : Total P", "Nitrate : Total P", "Inorganic NO3 : PO4"))
  )

# 3. Plot Heatmap (Now with 3 columns per Fjord Part)
p_heatmap <- ggplot(availability_comparison, aes(x = Year, y = Month, fill = Cluster_Count)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  # This facet_grid will now automatically generate the 3 ratio columns
  facet_grid(Fjord_Part ~ Ratio_Type) +
  theme_minimal() +
  labs(
    title = "Data Availability: Total vs. Hybrid vs. Inorganic Ratios",
    subtitle = "Comparing spatial coverage (500m clusters) of all three stoichiometry metrics",
    fill = "Samples", x = "Year", y = ""
  ) +
  theme(
    panel.grid = element_blank(), 
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_heatmap)

# ==============================================================================
# MODULE 3: SEASONAL STOICHIOMETRY TRENDS (All 3 Ratios)
# ==============================================================================

# 1. Aggregate and convert to Molarity
df_final_stoich <- oslo_final_sf %>%
  st_drop_geometry() %>%
  # ADDED: Phosphate included to calculate the purely inorganic ratio
  filter(parameter %in% c("Nitrogen", "Phosphorus", "Nitrate", "Phosphate")) %>%
  group_by(Fjord_Part, Year, Season, parameter) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    N_tot_uM = if("Nitrogen" %in% names(.)) Nitrogen / 14.01 else NA_real_,
    P_tot_uM = if("Phosphorus" %in% names(.)) Phosphorus / 30.97 else NA_real_,
    NO3_uM   = if("Nitrate" %in% names(.)) Nitrate / 14.01 else NA_real_,
    PO4_uM   = if("Phosphate" %in% names(.)) Phosphate / 30.97 else NA_real_,
    
    # Calculate all three ratios
    Ratio_TN_P    = N_tot_uM / P_tot_uM,
    Ratio_NO3_P   = NO3_uM / P_tot_uM,
    Ratio_NO3_PO4 = NO3_uM / PO4_uM
  ) %>%
  filter(!is.na(Season))

# 2. Format for Multi-Panel GAM Plot
df_trends_long <- df_final_stoich %>%
  select(Fjord_Part, Year, Season, Ratio_TN_P, Ratio_NO3_P, Ratio_NO3_PO4) %>%
  pivot_longer(cols = starts_with("Ratio"), names_to = "Ratio_Type", values_to = "Value") %>%
  # Cap at 100 to preserve Redfield visibility and remove NAs
  filter(is.finite(Value), Value >= 0, Value <= 100) %>% 
  mutate(Ratio_Type = recode(Ratio_Type, 
                             "Ratio_TN_P"    = "Total N : Total P", 
                             "Ratio_NO3_P"   = "Nitrate : Total P",
                             "Ratio_NO3_PO4" = "Inorganic NO3 : PO4")) %>%
  # Set order for the legend
  mutate(Ratio_Type = factor(Ratio_Type, levels = c("Total N : Total P", "Nitrate : Total P", "Inorganic NO3 : PO4")))

# 3. Plot Seasonal Trends
p_trends <- ggplot(df_trends_long, aes(x = Year, y = Value, color = Ratio_Type, fill = Ratio_Type)) +
  geom_hline(yintercept = 16, linetype = "dashed", color = "black", alpha = 0.8, linewidth = 0.8) +
  geom_point(alpha = 0.15, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1) +
  facet_grid(Season ~ Fjord_Part, scales = "free_y") +
  
  # UPDATED: 3 distinct colors (Green = Total, Blue = Hybrid, Purple = Inorganic)
  scale_color_manual(values = c("#4daf4a", "#377eb8", "#984ea3")) +
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#984ea3")) +
  
  theme_bw() +
  labs(
    title = "Oslofjord Nutrient Stoichiometry Trends (1930-2025)",
    subtitle = "Comparing Total Capacity vs. Immediate Bio-availability",
    y = "Molar Ratio (N:P)", x = "Year", color = "Ratio Type", fill = "Ratio Type"
  ) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

print(p_trends)

# 4. Extract Decadal Slopes
decadal_slopes <- df_trends_long %>%
  filter(!is.na(Value)) %>% # Clean out NAs first
  group_by(Fjord_Part, Season, Ratio_Type) %>%
  filter(n() >= 3) %>%      # THE FIX: Only run if there are >= 3 data points
  do(tidy(lm(Value ~ Year, data = .))) %>%
  filter(term == "Year") %>%
  mutate(Change_Per_Decade = estimate * 10) %>% 
  select(Fjord_Part, Season, Ratio_Type, estimate, p.value, Change_Per_Decade)

cat("\n[Module 3] Decadal Slopes Summary:\n")
print(decadal_slopes)


# ==============================================================================
# MODULE 5: BIOAVAILABLE FRACTIONS (NO3:TN and PO4:TP)
# Objective: Track what percentage of the Total pool is immediately bio-available.
# ==============================================================================

# 1. Aggregate, Convert, and Calculate Fractions
df_fractions <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(parameter %in% c("Nitrogen", "Phosphorus", "Nitrate", "Phosphate")) %>%
  
  # Group by Region, Year, and Season
  group_by(Fjord_Part, Year, Season, parameter) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  
  # Calculate molarity and fractions
  mutate(
    N_tot_uM = if("Nitrogen" %in% names(.)) Nitrogen / 14.01 else NA_real_,
    P_tot_uM = if("Phosphorus" %in% names(.)) Phosphorus / 30.97 else NA_real_,
    NO3_uM   = if("Nitrate" %in% names(.)) Nitrate / 14.01 else NA_real_,
    PO4_uM   = if("Phosphate" %in% names(.)) Phosphate / 30.97 else NA_real_,
    
    # Calculate the fractions (as percentages)
    Pct_NO3_TN = (NO3_uM / N_tot_uM),
    Pct_PO4_TP = (PO4_uM / P_tot_uM)
  ) %>%
  # Filter out rows with missing data or impossible fractions (e.g., > 100%)
  filter(!is.na(Season), Pct_NO3_TN <= 1.05, Pct_PO4_TP <= 1.05)

# 2. Format for Multi-Panel Plotting
df_frac_long <- df_fractions %>%
  select(Fjord_Part, Year, Season, Pct_NO3_TN, Pct_PO4_TP) %>%
  pivot_longer(cols = starts_with("Pct"), names_to = "Fraction_Type", values_to = "Percentage") %>%
  filter(is.finite(Percentage), Percentage >= 0) %>%
  mutate(
    Fraction_Type = recode(Fraction_Type, 
                           "Pct_NO3_TN" = "Nitrate (% of Total N)", 
                           "Pct_PO4_TP" = "Phosphate (% of Total P)"),
    # Set factor order
    Fraction_Type = factor(Fraction_Type, levels = c("Nitrate (% of Total N)", "Phosphate (% of Total P)"))
  )

# 3. Plotting the Fraction Trends
p_fractions <- ggplot(df_frac_long, aes(x = Year, y = Percentage, color = Fraction_Type, fill = Fraction_Type)) +
  
  # Raw data points (faded)
  geom_point(alpha = 0.2, size = 1) +
  
  # GAM Trendlines
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1.2) +
  
  # Facet by Season (Rows) and Fjord Part (Columns)
  facet_grid(Season ~ Fjord_Part) +
  
  # Convert Y-axis to neat percentages
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  
  # High-contrast colors for N vs P
  scale_color_manual(values = c("Nitrate (% of Total N)" = "#2166ac", "Phosphate (% of Total P)" = "#b2182b")) +
  scale_fill_manual(values = c("Nitrate (% of Total N)" = "#2166ac", "Phosphate (% of Total P)" = "#b2182b")) +
  
  theme_bw() +
  labs(
    title = "Ecosystem Quality: Shifting Bioavailable Fractions",
    subtitle = "Tracking the percentage of the Total nutrient pool that remains dissolved and inorganic",
    y = "Bioavailable Fraction (%)",
    x = "Year",
    color = "Nutrient Pool",
    fill = "Nutrient Pool"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

print(p_fractions)

# 4. Extract Decadal Slopes for the Text
fraction_slopes <- df_frac_long %>%
  filter(!is.na(Percentage)) %>% # Clean out NAs first
  group_by(Fjord_Part, Season, Fraction_Type) %>%
  filter(n() >= 3) %>%           # THE FIX
  do(tidy(lm(Percentage ~ Year, data = .))) %>%
  filter(term == "Year") %>%
  mutate(Change_Pct_Per_Decade = (estimate * 10) * 100) %>% 
  select(Fjord_Part, Season, Fraction_Type, estimate, p.value, Change_Pct_Per_Decade)

cat("\n[Module 5] Decadal Shifts in Bioavailable Fractions (% Change per Decade):\n")
print(fraction_slopes)

# ==============================================================================
# UPDATED MODULE 5: DEPTH-STRATIFIED BIOAVAILABLE FRACTIONS
# ==============================================================================

# 1. Aggregate and calculate fractions WITH Depth Bins
df_fractions_depth <- oslo_final_sf %>%
  st_drop_geometry() %>%
  mutate(
    Depth_Bin = case_when(
      depth <= 5 ~ "0-5m (Surface)",
      depth > 5 & depth <= 20 ~ "5-20m (Pycnocline/DCM)",
      TRUE ~ "Other" # We exclude deep water to focus on the photic zone
    )
  ) %>%
  filter(parameter %in% c("Nitrogen", "Phosphorus", "Nitrate", "Phosphate"), 
         Depth_Bin != "Other") %>%
  
  # Group by Depth_Bin as well!
  group_by(Fjord_Part, Depth_Bin, Year, Season, parameter) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  
  mutate(
    N_tot_uM = if("Nitrogen" %in% names(.)) Nitrogen / 14.01 else NA_real_,
    P_tot_uM = if("Phosphorus" %in% names(.)) Phosphorus / 30.97 else NA_real_,
    NO3_uM   = if("Nitrate" %in% names(.)) Nitrate / 14.01 else NA_real_,
    PO4_uM   = if("Phosphate" %in% names(.)) Phosphate / 30.97 else NA_real_,
    
    Pct_NO3_TN = (NO3_uM / N_tot_uM),
    Pct_PO4_TP = (PO4_uM / P_tot_uM)
  ) %>%
  filter(!is.na(Season), Pct_NO3_TN <= 1.05, Pct_PO4_TP <= 1.05)

# 2. Format for Plotting
df_frac_long <- df_fractions_depth %>%
  select(Fjord_Part, Depth_Bin, Year, Season, Pct_NO3_TN, Pct_PO4_TP) %>%
  pivot_longer(cols = starts_with("Pct"), names_to = "Fraction_Type", values_to = "Percentage") %>%
  mutate(
    Fraction_Type = recode(Fraction_Type, 
                           "Pct_NO3_TN" = "Nitrate (% of Total N)", 
                           "Pct_PO4_TP" = "Phosphate (% of Total P)"),
    Fraction_Type = factor(Fraction_Type, levels = c("Nitrate (% of Total N)", "Phosphate (% of Total P)"))
  )

# 3. Plotting: Now Faceted by Season and Depth
p_fractions_depth <- ggplot(df_frac_long, aes(x = Year, y = Percentage, color = Fraction_Type, fill = Fraction_Type)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1.2) +
  
  # THE CHANGE: We now facet by Season (rows) and Depth (columns)
  facet_grid(Season ~ Depth_Bin) + 
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Nitrate (% of Total N)" = "#2166ac", "Phosphate (% of Total P)" = "#b2182b")) +
  scale_fill_manual(values = c("Nitrate (% of Total N)" = "#2166ac", "Phosphate (% of Total P)" = "#b2182b")) +
  theme_bw() +
  labs(
    title = "Nutrient Quality in the Photic Zone",
    subtitle = "Comparing Bioavailable Fractions in the Surface vs. DCM",
    y = "Bioavailable Fraction (%)", x = "Year"
  ) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

print(p_fractions_depth)

# 1. Clean and Filter Data
df_frac_clean <- df_frac_long %>%
  # Remove Inner Oslofjord
  filter(Fjord_Part != "Inner Oslofjord") %>%
  # Ensure the factor order for regions is clean for the remaining three
  mutate(Fjord_Part = factor(Fjord_Part))

# 2. Define a Standard Plotting Function for consistency
plot_fraction_final <- function(data, depth_filter, plot_title) {
  ggplot(data %>% filter(Depth_Bin == depth_filter), 
         aes(x = Year, y = Percentage, color = Fraction_Type, fill = Fraction_Type)) +
    
    # Reference points (faded)
    geom_point(alpha = 0.15, size = 1.2) +
    
    # High-quality GAM Trendlines
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1.3) +
    
    # Facet Grid: Seasons (Rows) vs. Remaining 3 Regions (Columns)
    facet_grid(Season ~ Fjord_Part) + 
    
    # Styling
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = c("Nitrate (% of Total N)" = "#2166ac", "Phosphate (% of Total P)" = "#b2182b")) +
    scale_fill_manual(values = c("Nitrate (% of Total N)" = "#2166ac", "Phosphate (% of Total P)" = "#b2182b")) +
    
    theme_bw() +
    labs(
      title = plot_title,
      subtitle = paste("Tracking nutrient quality in the", depth_filter),
      y = "Bioavailable Fraction (%)",
      x = "Year",
      color = "Fraction Type",
      fill = "Fraction Type"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )
}

# 3. Generate the Two Individual Figures
p_frac_surface <- plot_fraction_final(df_frac_clean, "0-5m (Surface)", "Figure X: Surface Nutrient Quality")
p_frac_dcm     <- plot_fraction_final(df_frac_clean, "5-20m (Pycnocline/DCM)", "Figure Y: DCM Nutrient Quality")

# 4. View the plots
print(p_frac_surface)
print(p_frac_dcm)


# ==============================================================================
# Objective: Generate the depth-stratified stoichiometry ratios needed below!
# ==============================================================================

df_stoich_long <- oslo_final_sf %>%
  st_drop_geometry() %>%
  mutate(
    Depth_Bin = case_when(
      depth <= 5 ~ "0-5m (Surface)",
      depth > 5 & depth <= 20 ~ "5-20m (Pycnocline/DCM)",
      TRUE ~ "Other"
    )
  ) %>%
  filter(parameter %in% c("Nitrogen", "Phosphorus", "Nitrate", "Phosphate"), 
         Depth_Bin != "Other", !is.na(Season)) %>%
  group_by(Fjord_Part, Depth_Bin, Year, Season, parameter) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    N_tot_uM = if("Nitrogen" %in% names(.)) Nitrogen / 14.01 else NA_real_,
    P_tot_uM = if("Phosphorus" %in% names(.)) Phosphorus / 30.97 else NA_real_,
    NO3_uM   = if("Nitrate" %in% names(.)) Nitrate / 14.01 else NA_real_,
    PO4_uM   = if("Phosphate" %in% names(.)) Phosphate / 30.97 else NA_real_,
    
    Ratio_TN_P    = N_tot_uM / P_tot_uM,
    Ratio_NO3_P   = NO3_uM / P_tot_uM,
    Ratio_NO3_PO4 = NO3_uM / PO4_uM
  ) %>%
  select(Fjord_Part, Depth_Bin, Year, Season, Ratio_TN_P, Ratio_NO3_P, Ratio_NO3_PO4) %>%
  pivot_longer(cols = starts_with("Ratio"), names_to = "Ratio_Type", values_to = "Value") %>%
  filter(is.finite(Value), Value >= 0, Value <= 100) %>%
  mutate(
    Ratio_Type = recode(Ratio_Type, 
                        "Ratio_TN_P"    = "Total N : Total P", 
                        "Ratio_NO3_P"   = "Nitrate : Total P",
                        "Ratio_NO3_PO4" = "Inorganic NO3 : PO4"),
    Ratio_Type = factor(Ratio_Type, levels = c("Total N : Total P", "Nitrate : Total P", "Inorganic NO3 : PO4"))
  )


# ==============================================================================
# MODULE 6: SEASONAL PHENOLOGY IN THE PHOTIC ZONE
# Objective: Track stoichiometry dynamics in the Surface vs DCM across seasons.
# ==============================================================================

# 1. Filter for the biologically active layers and ensure factors are clean
df_seasonal_photic <- df_stoich_long %>%
  filter(
    Depth_Bin %in% c("0-5m (Surface)", "5-20m (Pycnocline/DCM)"),
    !is.na(Season)
  ) %>%
  # THE FIX: Use Norwegian season names to set the chronological top-to-bottom order
  mutate(Season = factor(Season, levels = c("Vinter", "Vår", "Sommer", "Høst")))

# 2. Plotting the Seasonal vs Depth Matrix
p_seasonal_depth <- ggplot(df_seasonal_photic, aes(x = Year, y = Value, color = Ratio_Type, fill = Ratio_Type)) +
  geom_hline(yintercept = 16, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_point(alpha = 0.15, size = 1, color = "gray50", shape = 16) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1.2) +
  facet_grid(Season ~ Depth_Bin, scales = "free_y") +
  scale_color_manual(values = c(
    "Total N : Total P" = "#4daf4a", 
    "Nitrate : Total P" = "#377eb8", 
    "Inorganic NO3 : PO4" = "#984ea3"
  )) +
  scale_fill_manual(values = c(
    "Total N : Total P" = "#4daf4a", 
    "Nitrate : Total P" = "#377eb8", 
    "Inorganic NO3 : PO4" = "#984ea3"
  )) +
  theme_bw() +
  labs(
    title = "System-Wide Seasonal Phenology of the Photic Zone",
    subtitle = "Tracking the biological consumption and winter reset of nutrients (1930-2025)",
    y = "Molar Ratio (N:P)",
    x = "Year",
    color = "Stoichiometric Metric",
    fill = "Stoichiometric Metric"
  ) +
  theme(
    legend.position = "bottom",
    # THE FIX: Explicitly call ggplot2::margin
    strip.text.x = element_text(face = "bold", size = 11, margin = ggplot2::margin(b = 5, t = 5)),
    strip.text.y = element_text(face = "bold", size = 11, margin = ggplot2::margin(l = 5, r = 5)),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

print(p_seasonal_depth)

# ==============================================================================
# MODULE 7: REGIONAL SEASONAL PHENOLOGY (Surface vs DCM)
# ==============================================================================

# 1. Base Filter (Active layers and clean seasons)
df_seasonal_photic <- df_stoich_long %>%
  filter(
    Depth_Bin %in% c("0-5m (Surface)", "5-20m (Pycnocline/DCM)"),
    !is.na(Season)
  ) %>%
  mutate(Season = factor(Season, levels = c("Vinter", "Vår", "Sommer", "Høst")))

# ---------------------------------------------------------
# FIGURE A: THE SURFACE LAYER (0-5m)
# ---------------------------------------------------------
df_surface <- df_seasonal_photic %>% filter(Depth_Bin == "0-5m (Surface)")

p_surface <- ggplot(df_surface, aes(x = Year, y = Value, color = Ratio_Type, fill = Ratio_Type)) +
  geom_hline(yintercept = 16, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_point(alpha = 0.15, size = 1, color = "gray50", shape = 16) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1.2) +
  facet_grid(Season ~ Fjord_Part, scales = "free_y") +
  scale_color_manual(values = c("Total N : Total P" = "#4daf4a", "Nitrate : Total P" = "#377eb8", "Inorganic NO3 : PO4" = "#984ea3")) +
  scale_fill_manual(values = c("Total N : Total P" = "#4daf4a", "Nitrate : Total P" = "#377eb8", "Inorganic NO3 : PO4" = "#984ea3")) +
  theme_bw() +
  labs(
    title = "Surface Layer (0-5m): Seasonal Nutrient Phenology",
    y = "Molar Ratio (N:P)", x = "Year", color = "Metric", fill = "Metric"
  ) +
  theme(
    legend.position = "bottom",
    # THE FIX: Explicitly call ggplot2::margin
    strip.text.x = element_text(face = "bold", size = 10, margin = ggplot2::margin(b = 4, t = 4)),
    strip.text.y = element_text(face = "bold", size = 10, margin = ggplot2::margin(l = 4, r = 4)),
    strip.background = element_rect(fill = "#e8f0fe")
  )

print(p_surface)

# ---------------------------------------------------------
# FIGURE B: THE PYCNOCLINE / DCM LAYER (5-20m)
# ---------------------------------------------------------
df_dcm <- df_seasonal_photic %>% filter(Depth_Bin == "5-20m (Pycnocline/DCM)")

p_dcm <- ggplot(df_dcm, aes(x = Year, y = Value, color = Ratio_Type, fill = Ratio_Type)) +
  geom_hline(yintercept = 16, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_point(alpha = 0.15, size = 1, color = "gray50", shape = 16) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1.2) +
  facet_grid(Season ~ Fjord_Part, scales = "free_y") +
  scale_color_manual(values = c("Total N : Total P" = "#4daf4a", "Nitrate : Total P" = "#377eb8", "Inorganic NO3 : PO4" = "#984ea3")) +
  scale_fill_manual(values = c("Total N : Total P" = "#4daf4a", "Nitrate : Total P" = "#377eb8", "Inorganic NO3 : PO4" = "#984ea3")) +
  theme_bw() +
  labs(
    title = "Pycnocline/DCM Layer (5-20m): Seasonal Nutrient Phenology",
    y = "Molar Ratio (N:P)", x = "Year", color = "Metric", fill = "Metric"
  ) +
  theme(
    legend.position = "bottom",
    # THE FIX: Explicitly call ggplot2::margin
    strip.text.x = element_text(face = "bold", size = 10, margin = ggplot2::margin(b = 4, t = 4)),
    strip.text.y = element_text(face = "bold", size = 10, margin = ggplot2::margin(l = 4, r = 4)),
    strip.background = element_rect(fill = "#e6f5d0") 
  )

print(p_dcm)

# ---------------------------------------------------------
# EXTRACTING THE REGIONAL SUMMER SLOPES
# ---------------------------------------------------------
summer_regional_slopes <- df_seasonal_photic %>%
  filter(Season == "Sommer", !is.na(Value)) %>%
  
  # THE FIX: We must group by Fjord_Part as well now!
  group_by(Fjord_Part, Depth_Bin, Ratio_Type) %>%
  
  filter(n() >= 3) %>% # Safety valve for empty historical buckets
  do(tidy(lm(Value ~ Year, data = .))) %>%
  filter(term == "Year") %>%
  mutate(Change_Per_Decade = estimate * 10) %>%
  select(Fjord_Part, Depth_Bin, Ratio_Type, estimate, p.value, Change_Per_Decade) %>%
  arrange(Depth_Bin, Fjord_Part)

cat("\n[Module 6/7] Critical Summer Shifts by Region and Depth:\n")
print(summer_regional_slopes, n = Inf)
