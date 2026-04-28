# ==============================================================================
# FIGURE 1 COMPILER: PHYSICAL FORCING (Layout Fixed)
# Assembles Salinity, Pycnocline (Depth & Strength), and Secchi into a composite A4 map
# ==============================================================================

library(tidyverse)
library(sf)
library(patchwork)

cat("\n--- 1. Generating Base Plots ---\n")

# A. Generate 'p_salinity' (Surface 0-5m)
p_salinity <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(grepl("Salin", parameter, ignore.case = TRUE), depth <= 5) %>%
  mutate(Year = year(date), value = as.numeric(value)) %>%
  filter(is.finite(value)) %>%
  group_by(Fjord_Part, Year) %>%
  summarise(Mean_Val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Year, y = Mean_Val)) +
  geom_point(alpha = 0.4, size = 1.5, color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              color = "#377eb8", linewidth = 1.2, alpha = 0.25) +
  facet_wrap(~ Fjord_Part, ncol = 4) +
  theme_bw() +
  labs(title = "Salinity Annual Mean (0-5m)", y = "Salinity (psu)", x = NULL) + # Dropped 'Year' to save space
  theme(strip.text = element_text(face = "bold", size = 10), panel.grid.minor = element_blank())

# B. Generate 'p_secchi' 
p_secchi <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(grepl("Secchi", parameter, ignore.case = TRUE)) %>%
  mutate(Year = year(date), value = as.numeric(value)) %>%
  filter(is.finite(value)) %>%
  group_by(Fjord_Part, Year) %>%
  summarise(Mean_Val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Year, y = Mean_Val)) +
  geom_point(alpha = 0.4, size = 1.5, color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              color = "#4daf4a", linewidth = 1.2, alpha = 0.25) +
  scale_y_reverse() + 
  facet_wrap(~ Fjord_Part, ncol = 4) +
  theme_bw() +
  labs(title = "Secchi Annual Mean", y = "Secchi Depth (m)", x = "Year") +
  theme(strip.text = element_text(face = "bold", size = 10), panel.grid.minor = element_blank())


cat("\n--- 2. Formatting Pycnocline Plots for Layout ---\n")

# Check to ensure your pycnocline boxplots exist in the environment
if(!exists("p_pyc_boxplot") | !exists("p_strength_boxplot")) {
  stop("ERROR: Pycnocline plots not found! Please run your pycnocline script first.")
}

# 1. Clean the Depth Boxplot: Strip subtitles, adjust margins, remove X-axis text so it locks onto the bottom plot
p_depth_clean <- p_pyc_boxplot +
  labs(title = "Pycnocline Depth Shift", subtitle = NULL, x = NULL, y = "Change in Depth\n(m / decade)") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    plot.margin = ggplot2::margin(b = 2) # Pulls it closer to the bottom plot
  )

# 2. Clean the Strength Boxplot: Strip subtitles, format Y-axis
p_strength_clean <- p_strength_boxplot +
  labs(title = "Pycnocline Strength Shift", subtitle = NULL, x = NULL, y = "Change in Strength\n(kg/m³ / m / decade)") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = ggplot2::margin(t = 2) # Pulls it closer to the top plot
  )

# 3. Stack them into Panel B, explicitly forcing the legend to the bottom
panel_B_pycnocline <- wrap_elements(
  (p_depth_clean / p_strength_clean) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom", legend.title = element_text(size = 9), legend.margin = ggplot2::margin(t = -10))
)


cat("\n--- 3. Compiling Master Figure 1 ---\n")

# Define layout (A=Salinity, C=Secchi, B=Pycnocline stack)
custom_layout <- "
  AB
  CB
"

# Assemble the master figure
figure_1 <- p_salinity + panel_B_pycnocline + p_secchi + 
  plot_layout(design = custom_layout, widths = c(1, 0.8)) + 
  plot_annotation(
    tag_levels = 'a', 
    tag_prefix = '(',
    tag_suffix = ')'
  ) & 
  theme(
    plot.tag = element_text(face = 'bold', size = 16),
    # THE FIX: vjust = 1 perfectly aligns the rotated text so it doesn't overlap the graph!
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

print(figure_1)

# Export to A4 Landscape
ggsave("Figure_1_Physical_Forcing_A4.pdf", figure_1, width = 11.69, height = 8.27, units = "in", dpi = 600)

cat("\nExport complete! Check your folder for the fixed 'Figure_1_Physical_Forcing_A4.pdf'\n")

# ==============================================================================
# BULLETPROOF FIGURE 2 COMPILER
# Fixes X-axis scaling in Panel C and adds 45-degree labels to Panel B
# ==============================================================================

library(tidyverse)
library(patchwork)

cat("\n--- 1. Regenerating TN & TP Panels ---\n")

# A. Prepare Data
df_nutrients_1990 <- df_analysis %>%
  filter(parameter %in% c("Total_Nitrogen", "Total_Phosphorus"), Year >= 1990) %>%
  group_by(parameter, Fjord_Part, Depth_Bin, Season) %>%
  mutate(
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE),
    IQR_val = Q3 - Q1,
    Lower_Bound = Q1 - 1.5 * IQR_val,
    Upper_Bound = Q3 + 1.5 * IQR_val
  ) %>%
  filter(value >= Lower_Bound | is.na(Lower_Bound), value <= Upper_Bound | is.na(Upper_Bound)) %>%
  ungroup()

annual_avg_nutrients <- df_nutrients_1990 %>%
  group_by(parameter, Fjord_Part, Depth_Bin, Year, Season) %>%
  summarise(Seasonal_Mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(parameter, Fjord_Part, Depth_Bin, Year) %>%
  summarise(Mean_Val = mean(Seasonal_Mean, na.rm = TRUE), n_seasons = n(), .groups = "drop")

# B. Plot Panel A (TN) - X-axis text remains blank to merge cleanly with Panel B
p_tn_1990 <- annual_avg_nutrients %>% 
  filter(parameter == "Total_Nitrogen") %>%
  ggplot(aes(x = Year, y = Mean_Val)) +
  geom_point(aes(size = n_seasons), color = "black", alpha = 0.6) +
  geom_smooth(method = "lm", color = "#e41a1c", linetype = "dashed", linewidth = 1, se = TRUE, alpha = 0.15) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=4), color = "#377eb8", linewidth = 1.2, se = FALSE) +
  scale_size_continuous(range = c(1.5, 4.5), breaks = 1:4, limits = c(1, 4)) +
  facet_grid(Depth_Bin ~ Fjord_Part, scales = "free_y") + 
  scale_x_continuous(breaks = seq(1990, 2025, by = 10), limits = c(1990, 2024)) +
  theme_bw() +
  labs(title = "Stratified Dynamics: Total Nitrogen (1990-Present)", y = "Total Nitrogen (µM)", x = NULL, size = "Seasons\nSampled") +
  theme(
    strip.text = element_text(face = "bold", size = 10), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# C. Plot Panel B (TP) - THE FIX: 45-degree angle applied to the Year ticks
p_tp_1990 <- annual_avg_nutrients %>% 
  filter(parameter == "Total_Phosphorus") %>%
  ggplot(aes(x = Year, y = Mean_Val)) +
  geom_point(aes(size = n_seasons), color = "black", alpha = 0.6) +
  geom_smooth(method = "lm", color = "#e41a1c", linetype = "dashed", linewidth = 1, se = TRUE, alpha = 0.15) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=4), color = "#377eb8", linewidth = 1.2, se = FALSE) +
  scale_size_continuous(range = c(1.5, 4.5), breaks = 1:4, limits = c(1, 4)) +
  facet_grid(Depth_Bin ~ Fjord_Part, scales = "free_y") + 
  scale_x_continuous(breaks = seq(1990, 2025, by = 10), limits = c(1990, 2024)) +
  theme_bw() +
  labs(title = "Stratified Dynamics: Total Phosphorus (1990-Present)", y = "Total Phosphorus (µM)", x = "Year") +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none",
    # THE FIX: Rotate Year text to 45 degrees
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) 
  )


cat("\n--- 2. Regenerating Stoichiometry Profile (Dynamic Scaling) ---\n")

# A. Prepare Data
ratio_slopes_all <- oslo_final_sf %>%
  st_drop_geometry() %>%
  mutate(Depth_Bin = case_when(
    depth <= 5 ~ "0-5m (Surface)",
    depth > 5 & depth <= 20 ~ "5-20m (Pycnocline/DCM)",
    depth > 20 & depth <= 50 ~ "20-50m",
    depth > 50 ~ ">50m (Deep Basin)",
    TRUE ~ "Unknown"
  )) %>%
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
  mutate(model = map(data, ~lm(Ratio_NO3_PO4 ~ Year, data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  filter(term == "Year") %>%
  mutate(
    Decadal_Ratio_Change = estimate * 10,
    Depth_Bin = factor(Depth_Bin, levels = c(">50m (Deep Basin)", "20-50m", "5-20m (Pycnocline/DCM)", "0-5m (Surface)"))
  ) %>%
  # THE FIX: Filter out the extreme non-physical anomalies (e.g., -230) so auto-scaling works!
  filter(Decadal_Ratio_Change > -30 & Decadal_Ratio_Change < 30)

# B. Plot Profile - Letting coord_flip() automatically encompass all valid data points
p_vertical_profile <- ggplot(ratio_slopes_all, aes(x = Depth_Bin, y = Decadal_Ratio_Change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  # THE FIX: Replaced hardcoded bounds with -Inf and Inf so backgrounds dynamically stretch
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#e8f0fe", alpha = 0.6) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#fce8e6", alpha = 0.6) +
  geom_boxplot(fill = "white", color = "black", alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_point(aes(fill = Fjord_Part), position = position_jitter(width = 0.15), size = 2.5, shape = 21, color = "black", alpha = 0.8) +
  
  # Removed the hardcoded ylim entirely!
  coord_flip() +
  
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "System-Wide Vertical Profile of Stoichiometric Shifts",
    subtitle = "Blue zone = Shifting toward N-Limitation. Red zone = Shifting toward P-Limitation.",
    y = "Change in NO3:PO4 Ratio (Units per Decade)", x = "Depth Horizon", fill = "Fjord Region"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom" 
  )


cat("\n--- 3. Compiling Master Figure 2 ---\n")

figure_2 <- p_tn_1990 / p_tp_1990 / p_vertical_profile + 
  plot_layout(heights = c(1, 1, 1.2)) + 
  plot_annotation(
    tag_levels = 'a', 
    tag_prefix = '(',
    tag_suffix = ')'
  ) & 
  theme(plot.tag = element_text(face = 'bold', size = 16))

print(figure_2)

# Export to A4 Portrait
ggsave("Figure_2_Nutrients_Stoichiometry_A4.pdf", figure_2, width = 8.27, height = 11.69, units = "in", dpi = 600)

cat("\nExport complete! Check your folder for the fully repaired 'Figure_2_Nutrients_Stoichiometry_A4.pdf'\n")

# ==============================================================================
# STRATIFIED SEASONAL STOICHIOMETRY: SURFACE VS. DCM
# Objective: Generate a faceted 4x4 plot comparing N:P ratios across both 
# depth layers to visualize vertical decoupling.
# ==============================================================================

library(tidyverse)

cat("\n--- Processing Stratified Seasonal Stoichiometry ---\n")

# 1. Extract and format data for both Surface and DCM layers
df_stoich_stratified <- df_analysis %>%
  filter(
    Depth_Bin %in% c("0-5m (Surface)", "5-20m (Pycnocline/DCM)"),
    parameter %in% c("Total_Nitrogen", "Total_Phosphorus", "Nitrate", "Phosphate")
  ) %>%
  group_by(Fjord_Part, Year, Season, Depth_Bin, parameter) %>%
  summarise(value_uM = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value_uM) %>%
  
  # 2. Calculate ratios
  mutate(
    `Total N : Total P` = Total_Nitrogen / Total_Phosphorus,
    `Nitrate : Total P` = Nitrate / Total_Phosphorus,
    `Inorganic NO3 : PO4` = Nitrate / Phosphate
  ) %>%
  
  # Filter extreme outliers for visualization clarity
  filter(
    if_all(contains(":"), ~ is.finite(.) & . < 150)
  ) %>%
  
  # Pivot to long format
  pivot_longer(
    cols = contains(":"),
    names_to = "Metric",
    values_to = "Ratio_Value"
  ) %>%
  
  # Ensure explicit ordering
  mutate(
    Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    Metric = factor(Metric, levels = c("Total N : Total P", "Nitrate : Total P", "Inorganic NO3 : PO4")),
    # Simplify label for legend
    Layer = ifelse(grepl("Surface", Depth_Bin), "Surface (0-5m)", "DCM (5-20m)")
  )

cat("\n--- Generating Stratified Seasonal Plot ---\n")

# 3. Create the Plot
p_stratified_stoich <- ggplot(df_stoich_stratified, aes(x = Year, y = Ratio_Value, color = Metric)) +
  
  # Redfield Ratio baseline
  geom_hline(yintercept = 16, linetype = "dotted", color = "black", linewidth = 0.8) +
  
  # Points with depth-specific shapes
  geom_point(aes(shape = Layer), alpha = 0.15, size = 1) +
  
  # GAM curves with linetype mapped to depth
  geom_smooth(aes(linetype = Layer, fill = Metric), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              linewidth = 1, alpha = 0.1) +
  
  facet_grid(Season ~ Fjord_Part) +
  
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_x_continuous(breaks = seq(1980, 2020, by = 20)) +
  
  scale_color_manual(values = c("Total N : Total P" = "#4daf4a", 
                                "Nitrate : Total P" = "#984ea3", 
                                "Inorganic NO3 : PO4" = "#377eb8")) +
  scale_fill_manual(values = c("Total N : Total P" = "#4daf4a", 
                               "Nitrate : Total P" = "#984ea3", 
                               "Inorganic NO3 : PO4" = "#377eb8")) +
  
  theme_bw() +
  labs(
    title = "Stratified Seasonal Nutrient Stoichiometry",
    subtitle = "Solid lines = Surface (0-5m), Dashed lines = DCM (5-20m)",
    y = "Molar Ratio (N:P)",
    x = "Year",
    linetype = "Depth Horizon",
    shape = "Depth Horizon"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_stratified_stoich)

# Export
ggsave("Supplementary_Stratified_Stoich_Seasonal_A4.pdf", p_stratified_stoich, width = 11, height = 9)

# ==============================================================================
# STRATIFIED SEASONAL NUTRIENTS (DUAL Y-AXIS - 10X MAGNIFICATION)
# Scaling: Left Axis (N) = 0-60 µM | Right Axis (P) = 0-6 µM
# ==============================================================================

library(tidyverse)

# 1. Prepare and Transform Data
df_nutrients_6_60 <- df_analysis %>%
  filter(
    Depth_Bin %in% c("0-5m (Surface)", "5-20m (Pycnocline/DCM)"),
    parameter %in% c("Total_Nitrogen", "Total_Phosphorus", "Nitrate", "Phosphate"),
    !is.na(Season), !is.na(Fjord_Part)
  ) %>%
  # Normalize Season names to Norwegian labels
  mutate(Season = case_when(
    Season %in% c("Winter", "Vinter") ~ "Vinter",
    Season %in% c("Spring", "Vår")    ~ "Vår",
    Season %in% c("Summer", "Sommer") ~ "Sommer",
    Season %in% c("Autumn", "Høst")   ~ "Høst",
    TRUE ~ Season
  )) %>%
  mutate(
    Season = factor(Season, levels = c("Vinter", "Vår", "Sommer", "Høst")),
    Layer = ifelse(grepl("Surface", Depth_Bin), "Surface", "DCM")
  ) %>%
  group_by(Fjord_Part, Year, Season, Layer, parameter) %>%
  summarise(value_uM = mean(value, na.rm = TRUE), .groups = "drop") %>%
  # Transform Phosphorus values by 10x to align 6 uM with the 60 uM N-limit
  mutate(plot_value = ifelse(parameter %in% c("Total_Phosphorus", "Phosphate"), 
                             value_uM * 10, 
                             value_uM))

# 2. Generate the 16-Panel Plot
p_nutrients_6_60 <- ggplot(df_nutrients_6_60, aes(x = Year, y = plot_value, color = parameter)) +
  geom_point(aes(shape = Layer), alpha = 0.15, size = 0.7) +
  geom_smooth(aes(linetype = Layer, fill = parameter), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              linewidth = 0.8, alpha = 0.15) +
  
  facet_grid(Season ~ Fjord_Part) +
  
  # Dual Axis Configuration: 0-60 Primary, 0-6 Secondary
  scale_y_continuous(
    name = "Nitrogen (TN / Nitrate) [µM]",
    limits = c(0, 60),
    breaks = seq(0, 60, by = 10),
    sec.axis = sec_axis(~ . / 10, name = "Phosphorus (TP / Phosphate) [µM]")
  ) +
  
  scale_color_manual(values = c("Total_Nitrogen" = "#e41a1c", "Nitrate" = "#377eb8", 
                                "Total_Phosphorus" = "#4daf4a", "Phosphate" = "#984ea3")) +
  scale_fill_manual(values = c("Total_Nitrogen" = "#e41a1c", "Nitrate" = "#377eb8", 
                               "Total_Phosphorus" = "#4daf4a", "Phosphate" = "#984ea3")) +
  theme_bw() +
  labs(
    title = "Surface vs. DCM Nutrient Dynamics: 10x Stoichiometric Scaling",
    subtitle = "Nitrogen scaled 0-60 µM | Phosphorus scaled 0-6 µM (10:1 visual ratio)",
    x = "Year"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    axis.title.y.right = element_text(color = "#4daf4a", face = "bold"),
    axis.title.y.left = element_text(color = "#e41a1c", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

# 3. Export to PDF
ggsave("Supplementary_16_Panel_Nutrients_6_60_Scale.pdf", 
       plot = p_nutrients_6_60, width = 14, height = 10)

cat("File successfully saved: Supplementary_16_Panel_Nutrients_6_60_Scale.pdf\n")

# Install patchwork if you haven't already
# install.packages("patchwork")
library(patchwork)

######## DATA AVAILABILITY PLOTS ################

# --- 1. Assign your first plot to a variable (p1) ---
p1 <- ggplot(coverage_check, aes(x = Decade, y = Fjord_Part, fill = n)) +
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

# --- 2. Assign your second plot to a variable (p2) ---
p2 <- oslo_final_sf %>%
  filter(depth <= 400) %>% 
  ggplot(aes(x = depth, fill = parameter)) +
  geom_histogram(binwidth = 5, alpha = 0.8, show.legend = FALSE) +
  scale_y_log10() +
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

# --- 3. Combine them using patchwork ---
# The "/" operator stacks them vertically. 
# You can use "|" to place them side-by-side if preferred.
combined_plot <- p1 / p2 + 
  plot_annotation(tag_levels = 'A') & # This automatically adds 'A' and 'B' tags
  theme(plot.tag = element_text(size = 18, face = 'bold'))

# --- 4. Save the combined high-resolution figure as a PDF ---
ggsave(
  filename = "Supplementary_fig_2_combined.pdf", 
  plot = combined_plot, 
  device = "pdf",
  width = 12,    # Adjust width as needed for your facet_wraps
  height = 14,   # Make it tall enough to fit both plots comfortably
  bg = "white"   # Ensures the background isn't transparent
)

