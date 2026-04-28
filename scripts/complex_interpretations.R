# ==============================================================================
# OSLOFJORD GRAND MASTER SCRIPT: PHYSICS, BIOGEOCHEMISTRY & COUPLING
# ==============================================================================

# --- Load Required Libraries ---
library(tidyverse)
library(ggplot2)
library(broom)
library(oce)       # For Seawater Density (Sigma-t) calculations
library(patchwork)
library(sf)

# ==============================================================================
# MODULE 1: TRUE PYCNOCLINE STRATIFICATION (PHYSICS)
# ==============================================================================

# 1A. Calculate Density and Extract the Maximum Gradient (Pycnocline)
df_true_pycnocline <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(parameter %in% c("Salinity", "Temperature"), !is.na(value), !is.na(depth)) %>%
  group_by(Fjord_Part, Year, Season, Week_Date, site_cluster, depth, parameter) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  filter(!is.na(Salinity), !is.na(Temperature)) %>%
  mutate(Density = swSigmaT(salinity = Salinity, temperature = Temperature, pressure = depth)) %>%
  group_by(Fjord_Part, Year, Season, Week_Date, site_cluster) %>%
  arrange(depth, .by_group = TRUE) %>%
  filter(n_distinct(depth) >= 2) %>%
  mutate(
    delta_Rho = Density - lag(Density),
    delta_Z = depth - lag(depth),
    Density_Gradient = delta_Rho / delta_Z
  ) %>%
  filter(!is.na(Density_Gradient), is.finite(Density_Gradient), delta_Z > 0) %>%
  slice_max(Density_Gradient, n = 1, with_ties = FALSE) %>%
  select(Fjord_Part, Year, Season, Week_Date, site_cluster, 
         Pycnocline_Depth = depth, Pycnocline_Strength = Density_Gradient) %>%
  ungroup()

# 1B. Extract Pycnocline Depth Slopes
pyc_slopes <- df_true_pycnocline %>%
  group_by(site_cluster) %>%
  filter(n_distinct(Year) >= 5) %>% 
  group_by(Fjord_Part, site_cluster) %>%
  nest() %>% 
  mutate(
    model = map(data, ~lm(Pycnocline_Depth ~ Year, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Year") %>%
  mutate(
    Decadal_Change = estimate * 10,
    Significant = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not Significant")
  )

# ==============================================================================
# MODULE 2: ABSOLUTE NUTRIENT TRAJECTORIES (RAW µM)
# ==============================================================================

# 2A. Prepare Raw Nutrient Data for the DCM
df_raw_dcm <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(depth > 5 & depth <= 20, parameter %in% c("Nitrate", "Phosphate")) %>%
  mutate(
    value_uM = case_when(
      parameter == "Nitrate" ~ as.numeric(value) / 14.01,
      parameter == "Phosphate" ~ as.numeric(value) / 30.97,
      TRUE ~ as.numeric(value)
    )
  ) %>%
  filter(is.finite(value_uM), value_uM > 0) %>%
  group_by(Fjord_Part, site_cluster, parameter, Year) %>%
  summarise(value_uM = mean(value_uM, na.rm = TRUE), .groups = "drop")

# 2B. Calculate Absolute Slopes
raw_nutrient_slopes <- df_raw_dcm %>%
  group_by(Fjord_Part, site_cluster, parameter) %>%
  filter(n_distinct(Year) >= 5) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(value_uM ~ Year, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Year") %>%
  mutate(Decadal_Change_uM = estimate * 10) %>%
  select(Fjord_Part, site_cluster, parameter, Decadal_Change_uM) %>%
  pivot_wider(names_from = parameter, values_from = Decadal_Change_uM) %>%
  rename(Slope_NO3_uM = Nitrate, Slope_PO4_uM = Phosphate) %>%
  filter(!is.na(Slope_NO3_uM), !is.na(Slope_PO4_uM))

# ==============================================================================
# MODULE 3: STOICHIOMETRIC RATIO TRAJECTORIES IN THE DCM
# ==============================================================================

# 3A. Prepare Site-Level Ratio Data
df_site_trends_dcm <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(depth > 5 & depth <= 20, parameter %in% c("Nitrate", "Phosphate")) %>%
  group_by(Fjord_Part, site_cluster, Year, date, parameter) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    NO3_uM = if("Nitrate" %in% names(.)) Nitrate / 14.01 else NA_real_,
    PO4_uM = if("Phosphate" %in% names(.)) Phosphate / 30.97 else NA_real_,
    Ratio_NO3_PO4 = NO3_uM / PO4_uM
  ) %>%
  filter(is.finite(Ratio_NO3_PO4), Ratio_NO3_PO4 <= 100) %>%
  group_by(site_cluster) %>%
  filter(n_distinct(Year) >= 5) %>%
  ungroup()

# 3B. Calculate DCM Ratio Slopes
dcm_slopes <- df_site_trends_dcm %>%
  group_by(Fjord_Part, site_cluster) %>%
  nest() %>% 
  mutate(
    model = map(data, ~lm(Ratio_NO3_PO4 ~ Year, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Year") %>%
  mutate(
    Decadal_Change = estimate * 10,
    Significant = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not Significant")
  )

# ==============================================================================
# MODULE 4: PHYSICAL-BIOGEOCHEMICAL COUPLING
# ==============================================================================

# 4A. Merge Pycnocline and DCM Ratio Slopes
coupled_slopes <- inner_join(
  dcm_slopes %>% select(Fjord_Part, site_cluster, Slope_NO3 = Decadal_Change, Sig_NO3 = Significant),
  pyc_slopes %>% select(Fjord_Part, site_cluster, Slope_Pyc = Decadal_Change, Sig_Pyc = Significant),
  by = c("Fjord_Part", "site_cluster")
) %>%
  mutate(
    Both_Significant = ifelse(grepl("p < 0.05", Sig_NO3) & grepl("p < 0.05", Sig_Pyc), 
                              "Both Significant", "One or Neither")
  )

# 4B. Plot the Coupling Scatterplot
p_coupled <- ggplot(coupled_slopes, aes(x = Slope_Pyc, y = Slope_NO3)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(fill = Fjord_Part, shape = Both_Significant), size = 3, alpha = 0.8) +
  scale_shape_manual(values = c("Both Significant" = 21, "One or Neither" = 24)) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE, alpha = 0.2) +
  scale_x_reverse() +
  facet_wrap(~Fjord_Part, scales = "free") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(
    title = "Coupling of Physical and Biogeochemical Shifts",
    subtitle = "Top-Right: Pycnocline getting shallower AND NO3:PO4 ratio increasing.",
    x = "Rate of Pycnocline Shift (M/Decade) -> [Negative = Shallower]",
    y = "Rate of NO3:PO4 Shift (Units/Decade) -> [Positive = Increasing]",
    fill = "Region", shape = "Statistical Rigor"
  ) +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

# ==============================================================================
# MODULE 5: VERTICAL PROFILES OF STOICHIOMETRY (ALL DEPTHS)
# ==============================================================================

# 5A. Calculate Ratio Slopes for Every Depth Horizon
ratio_slopes_all <- oslo_final_sf %>%
  st_drop_geometry() %>%
  mutate(
    Depth_Bin = case_when(
      depth <= 5 ~ "0-5m (Surface)",
      depth > 5 & depth <= 20 ~ "5-20m (Pycnocline/DCM)",
      depth > 20 & depth <= 50 ~ "20-50m",
      depth > 50 ~ ">50m (Deep Basin)",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(parameter %in% c("Nitrate", "Phosphate"), Depth_Bin != "Unknown") %>%
  group_by(Fjord_Part, site_cluster, Depth_Bin, parameter, Year) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    NO3_uM = if("Nitrate" %in% names(.)) Nitrate / 14.01 else NA_real_,
    PO4_uM = if("Phosphate" %in% names(.)) Phosphate / 30.97 else NA_real_,
    Ratio_NO3_PO4 = NO3_uM / PO4_uM
  ) %>%
  filter(is.finite(Ratio_NO3_PO4), Ratio_NO3_PO4 <= 100) %>%
  group_by(Fjord_Part, site_cluster, Depth_Bin) %>%
  filter(n_distinct(Year) >= 5) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(Ratio_NO3_PO4 ~ Year, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Year") %>%
  mutate(
    Decadal_Ratio_Change = estimate * 10,
    Significant = p.value < 0.05,
    Depth_Bin = factor(Depth_Bin, levels = c(">50m (Deep Basin)", "20-50m", "5-20m (Pycnocline/DCM)", "0-5m (Surface)"))
  )

# 5B. Plot System-Wide Vertical Consensus
p_vertical_profile <- ggplot(ratio_slopes_all, aes(x = Depth_Bin, y = Decadal_Ratio_Change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#e8f0fe", alpha = 0.6) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#fce8e6", alpha = 0.6) +
  geom_boxplot(fill = "white", color = "black", alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_point(aes(fill = Fjord_Part), position = position_jitter(width = 0.15), 
             size = 2.5, shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "System-Wide Vertical Profile of Stoichiometric Shifts",
    subtitle = "Blue zone = Shifting toward N-Limitation. Red zone = Shifting toward P-Limitation.",
    y = "Change in NO3:PO4 Ratio (Units per Decade)", x = "Depth Horizon", fill = "Fjord Region"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

# ==============================================================================
# PRINT ALL FIGURES FOR REVIEW
# ==============================================================================
print(p_coupled)
print(p_vertical_profile)

# To view the Redfield / Option plots from earlier, you can also run:
# print(p_redfield_global) 
# print(p_dcm_slopes)


