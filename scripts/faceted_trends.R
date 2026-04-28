# ==============================================================================
# MODULE 14: ENHANCED LONG-TERM TREND ANALYSIS
# Objective: Extract multi-decadal trends for all parameters (including TN/TP),
# split by Surface vs. DCM, with strict seasonal-bias removal.
# ==============================================================================

library(tidyverse)
library(sf)
library(broom)

# ---------------------------------------------------------
# 1. Preparation (Bulletproof Name Matching)
# ---------------------------------------------------------
df_analysis <- oslo_final_sf %>%
  st_drop_geometry() %>%
  
  # Standardize all parameter names FIRST using our aggressive regex
  mutate(
    parameter_clean = case_when(
      grepl("Chl", parameter, ignore.case = TRUE) ~ "Chlorophyll-a",
      grepl("Nitrate", parameter, ignore.case = TRUE) ~ "Nitrate",
      grepl("Phosphate", parameter, ignore.case = TRUE) ~ "Phosphate",
      grepl("Tot.*N|Nitrogen", parameter, ignore.case = TRUE) ~ "Total_Nitrogen", 
      grepl("Tot.*P|Phosphorus", parameter, ignore.case = TRUE) ~ "Total_Phosphorus", 
      grepl("Temp", parameter, ignore.case = TRUE) ~ "Temperature",
      grepl("Salin", parameter, ignore.case = TRUE) ~ "Salinity",
      grepl("Secchi", parameter, ignore.case = TRUE) ~ "Secchi",
      TRUE ~ "IGNORE"
    )
  ) %>%
  
  # Filter out unwanted parameters
  filter(parameter_clean != "IGNORE") %>%
  
  mutate(
    Year = year(date),
    Depth_Bin = case_when(
      depth <= 5 ~ "0-5m (Surface)",
      depth > 5 & depth <= 20 ~ "5-20m (Pycnocline/DCM)",
      TRUE ~ "Deep"
    ),
    # Convert all Nitrogen/Phosphorus pools to µM
    value_clean = case_when(
      parameter_clean %in% c("Nitrate", "Total_Nitrogen") ~ as.numeric(value) / 14.01,
      parameter_clean %in% c("Phosphate", "Total_Phosphorus") ~ as.numeric(value) / 30.97,
      TRUE ~ as.numeric(value)
    )
  ) %>%
  
  # Drop deep data, keep valid numbers, ensure we have Season data
  filter(!is.na(Fjord_Part), !is.na(Season), Depth_Bin != "Deep", is.finite(value_clean), value_clean > 0) %>%
  
  # OVERRIDE the old 'parameter' and 'value' columns so the loop works perfectly
  mutate(
    parameter = parameter_clean,
    value = value_clean
  )

# Define the master list of parameters to loop through
params_to_analyze <- c("Temperature", "Salinity", "Secchi", "Nitrate", "Phosphate", "Total_Nitrogen", "Total_Phosphorus", "Chlorophyll-a")

# ---------------------------------------------------------
# 2. The Master Analysis Function
# ---------------------------------------------------------
analyze_depth_regional_trends <- function(p_name) {
  
  cat("\nProcessing:", p_name, "...\n")
  
  # 1. Filter by Parameter
  p_df <- df_analysis %>% filter(parameter == p_name)
  if(nrow(p_df) == 0) return(NULL) 
  
  # 2. Strict Outlier Removal (Protects true spring blooms!)
  p_df_clean <- p_df %>%
    group_by(Fjord_Part, Depth_Bin, Season) %>%
    mutate(
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      IQR_val = Q3 - Q1,
      Lower_Bound = Q1 - 1.5 * IQR_val,
      Upper_Bound = Q3 + 1.5 * IQR_val
    ) %>%
    filter(value >= Lower_Bound | is.na(Lower_Bound), 
           value <= Upper_Bound | is.na(Upper_Bound)) %>%
    ungroup()
  
  # 3. Bias-Free Averaging (Season -> Year) WITH DENSITY COUNT
  seasonal_avg <- p_df_clean %>%
    group_by(Fjord_Part, Depth_Bin, Year, Season) %>%
    summarise(Seasonal_Mean = mean(value, na.rm = TRUE), .groups = "drop")
  
  annual_avg <- seasonal_avg %>%
    group_by(Fjord_Part, Depth_Bin, Year) %>%
    summarise(
      Mean_Val = mean(Seasonal_Mean, na.rm = TRUE),
      n_seasons = n(), # Count how many distinct seasons made up this year
      .groups = "drop"
    )
  
  # 4. Generate the Faceted Plot
  plot <- ggplot(annual_avg, aes(x = Year, y = Mean_Val)) +
    
    # Yearly mean points (Size mapped to n_seasons)
    geom_point(aes(size = n_seasons), color = "black", alpha = 0.6) +
    
    # Linear Trend (Red Dashed)
    geom_smooth(method = "lm", color = "#e41a1c", linetype = "dashed", linewidth = 1, se = TRUE, alpha = 0.15) +
    
    # Non-linear Trend (Blue GAM)
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=4), 
                color = "#377eb8", linewidth = 1.2, se = FALSE) +
    
    scale_size_continuous(range = c(1.5, 4.5), breaks = 1:4, limits = c(1, 4)) +
    
    # Split by Depth (Rows) and Region (Columns)
    facet_grid(Depth_Bin ~ Fjord_Part, scales = "free_y") + 
    
    scale_x_continuous(breaks = seq(1930, 2025, by = 20)) +
    theme_bw() +
    labs(
      title = paste("Long-Term Stratified Dynamics:", p_name),
      subtitle = "Point size reflects data density (1 to 4 seasons sampled per year)",
      y = paste(p_name, "Annual Mean"),
      x = "Year",
      size = "Seasons\nSampled"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )
  
  print(plot)
  
  # 5. Extract Statistics
  stats <- annual_avg %>%
    group_by(Fjord_Part, Depth_Bin) %>%
    filter(n_distinct(Year) >= 10) %>%
    do(tidy(lm(Mean_Val ~ Year, data = .))) %>%
    filter(term == "Year") %>%
    mutate(
      parameter = p_name,
      Decadal_Shift = estimate * 10,
      Significant = ifelse(p.value < 0.05, "Yes", "No")
    ) %>%
    select(parameter, Fjord_Part, Depth_Bin, Decadal_Shift, p_value = p.value, Significant)
  
  return(stats)
}

# ---------------------------------------------------------
# 3. Execution and Export
# ---------------------------------------------------------
# Run the loop for all 8 parameters!
final_stats_list <- lapply(params_to_analyze, analyze_depth_regional_trends)

# Bind the results into a single master table
master_summary_table <- bind_rows(final_stats_list)

# Print the statistically significant findings to the console
cat("\n========================================================\n")
cat(" MASTER STATISTICAL SUMMARY (Significant Trends Only)\n")
cat("========================================================\n")
print(master_summary_table %>% filter(Significant == "Yes") %>% arrange(parameter, Depth_Bin), n = 50)

# Export to a CSV so you can format it nicely in Word or Excel
write_csv(master_summary_table, "Oslofjord_Master_Trends_Stratified.csv")
cat("\nData successfully exported to 'Oslofjord_Master_Trends_Stratified.csv'\n")

library(tidyverse)
library(patchwork)

cat("\n--- 2. Plotting Dual-Scale Nutrient Trends ---\n")

# A. Reshape the data so N and P are in separate columns for dual plotting
df_wide <- df_nutrients_abs %>%
  pivot_wider(names_from = parameter_clean, values_from = value_uM)

# Define the scaling factor to bring P up to the visual level of N
# 16 is the standard Redfield ratio, but you can adjust this if needed
scale_factor <- 16 

# ---------------------------------------------------------
# B. Bioavailable Nutrients (Dual Axis)
# ---------------------------------------------------------
p_bio_dual <- df_wide %>%
  drop_na(Nitrate, Phosphate) %>%
  ggplot(aes(x = Year)) +
  
  # Nitrate (Maps to Primary Left Y-Axis)
  geom_smooth(aes(y = Nitrate, color = "Nitrate", linetype = Depth_Bin), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              se = TRUE, alpha = 0.15, linewidth = 1.2) +
  
  # Phosphate (Multiplied by scale_factor to lift it visually)
  geom_smooth(aes(y = Phosphate * scale_factor, color = "Phosphate", linetype = Depth_Bin), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              se = TRUE, alpha = 0.15, linewidth = 1.2) +
  
  # Set up the Dual Y-Axis
  scale_y_continuous(
    name = "Nitrate (µM)", 
    sec.axis = sec_axis(~ . / scale_factor, name = "Phosphate (µM)")
  ) +
  
  facet_grid(Season ~ Fjord_Part, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("Nitrate" = "#2166ac", "Phosphate" = "#b2182b")) +
  scale_linetype_manual(values = c("0-5m (Surface)" = "solid", "5-20m" = "11")) +
  labs(
    title = "Absolute Bioavailable Nutrients (Dual Scale)",
    subtitle = "Solid lines = Surface (0-5m), Dotted lines = DCM (5-20m)",
    x = "Year", color = "Nutrient", linetype = "Depth Layer"
  ) +
  theme(
    strip.text = element_text(face = "bold"), 
    legend.position = "bottom",
    # Color the axis titles to match the lines for easy reading
    axis.title.y.left = element_text(color = "#2166ac", face = "bold"),
    axis.title.y.right = element_text(color = "#b2182b", face = "bold")
  )

print(p_bio_dual)

# ---------------------------------------------------------
# C. Total Nutrients (Dual Axis)
# ---------------------------------------------------------
p_total_dual <- df_wide %>%
  drop_na(Total_Nitrogen, Total_Phosphorus) %>%
  ggplot(aes(x = Year)) +
  
  # Total Nitrogen (Maps to Primary Left Y-Axis)
  geom_smooth(aes(y = Total_Nitrogen, color = "Total N", linetype = Depth_Bin), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              se = TRUE, alpha = 0.15, linewidth = 1.2) +
  
  # Total Phosphorus (Multiplied by scale_factor)
  geom_smooth(aes(y = Total_Phosphorus * scale_factor, color = "Total P", linetype = Depth_Bin), 
              method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              se = TRUE, alpha = 0.15, linewidth = 1.2) +
  
  # Set up the Dual Y-Axis
  scale_y_continuous(
    name = "Total Nitrogen (µM)", 
    sec.axis = sec_axis(~ . / scale_factor, name = "Total Phosphorus (µM)")
  ) +
  
  facet_grid(Season ~ Fjord_Part, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("Total N" = "#053061", "Total P" = "#67001f")) +
  scale_linetype_manual(values = c("0-5m (Surface)" = "solid", "5-20m" = "11")) +
  labs(
    title = "Absolute Total Nutrients (Dual Scale)",
    subtitle = "Solid lines = Surface (0-5m), Dotted lines = DCM (5-20m)",
    x = "Year", color = "Nutrient", linetype = "Depth Layer"
  ) +
  theme(
    strip.text = element_text(face = "bold"), 
    legend.position = "bottom",
    # Color the axis titles to match the lines for easy reading
    axis.title.y.left = element_text(color = "#053061", face = "bold"),
    axis.title.y.right = element_text(color = "#67001f", face = "bold")
  )

print(p_total_dual)

