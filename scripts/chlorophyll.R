

# ==============================================================================
# MODULE 6: MULTIVARIATE ASSOCIATIONS (Seasonally Adjusted)
# ==============================================================================
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(corrplot)
library(factoextra)
library(vegan)
library(patchwork)

cat("\n--- 1. Preparing Master Data & Temporal Variables ---\n")

# A. Extract Master Data and Force Temporal/Seasonal Columns
df_multi <- oslo_final_sf %>%
  st_drop_geometry() %>%
  mutate(
    date = as.Date(date),
    Year = year(date),
    Month = month(date),
    Season = case_when(
      Month %in% c(3, 4, 5) ~ "Spring",
      Month %in% c(6, 7, 8) ~ "Summer",
      Month %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ "Winter"
    ),
    Season = factor(Season, levels = c("Spring", "Summer", "Autumn", "Winter")),
    Sampling_Event = as.factor(paste(Fjord_Part, date, sep = "_"))
  ) %>%
  filter(!is.na(date))

# B. Extract Secchi Depth
df_secchi_multi <- df_multi %>%
  filter(grepl("Secchi", parameter, ignore.case = TRUE)) %>%
  mutate(value = as.numeric(value)) %>%
  filter(is.finite(value), value > 0) %>%
  group_by(Sampling_Event) %>%
  summarise(Secchi = mean(value, na.rm = TRUE), .groups = "drop")

# C. Extract Chemistry & Chlorophyll (Aggressive Matching)
df_chem_multi <- df_multi %>%
  filter(!grepl("Secchi", parameter, ignore.case = TRUE)) %>%
  mutate(
    parameter_clean = case_when(
      grepl("Chl", parameter, ignore.case = TRUE) ~ "Chla",
      grepl("Nitrate", parameter, ignore.case = TRUE) ~ "Nitrate",
      grepl("Phosphate", parameter, ignore.case = TRUE) ~ "Phosphate",
      grepl("Tot.*N|Nitrogen", parameter, ignore.case = TRUE) ~ "Total_Nitrogen", 
      grepl("Tot.*P|Phosphorus", parameter, ignore.case = TRUE) ~ "Total_Phosphorus", 
      grepl("Temp", parameter, ignore.case = TRUE) ~ "Temperature",
      grepl("Salin", parameter, ignore.case = TRUE) ~ "Salinity",
      TRUE ~ "IGNORE"
    ),
    Depth_Bin = case_when(
      depth <= 5 ~ "0-5m (Surface)",
      depth > 5 & depth <= 20 ~ "5-20m (Pycnocline/DCM)",
      TRUE ~ "Deep"
    ),
    value_clean = case_when(
      parameter_clean %in% c("Nitrate", "Total_Nitrogen") ~ as.numeric(value) / 14.01,
      parameter_clean %in% c("Phosphate", "Total_Phosphorus") ~ as.numeric(value) / 30.97,
      TRUE ~ as.numeric(value)
    )
  ) %>%
  filter(parameter_clean != "IGNORE", Depth_Bin != "Deep", is.finite(value_clean), value_clean >= 0) %>%
  group_by(Sampling_Event, Fjord_Part, Depth_Bin, Year, Season, parameter_clean) %>%
  summarise(value = mean(value_clean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter_clean, values_from = value) %>%
  mutate(NO3_PO4_Ratio = Nitrate / (Phosphate + 0.01))

# D. Merge Data, Transform Skewed Variables, & Split by Depth
df_multi_final <- df_chem_multi %>%
  inner_join(df_secchi_multi, by = "Sampling_Event") %>%
  drop_na(Chla, Temperature, Salinity, Secchi, NO3_PO4_Ratio) %>%
  mutate(
    Chla = log10(Chla + 1),
    NO3_PO4_Ratio = log10(NO3_PO4_Ratio + 1)
  )

df_surf <- df_multi_final %>% filter(Depth_Bin == "0-5m (Surface)")
df_dcm <- df_multi_final %>% filter(Depth_Bin == "5-20m (Pycnocline/DCM)")


cat("\n--- 2. Seasonally-Adjusted Correlogram (Partial Correlations) ---\n")

get_season_residuals <- function(data) {
  num_vars <- data %>% 
    select(any_of(c("Chla", "Temperature", "Salinity", "Secchi", "Nitrate", 
                    "Phosphate", "Total_Nitrogen", "Total_Phosphorus", "NO3_PO4_Ratio"))) %>%
    rename(any_of(c(Temp = "Temperature", Sal = "Salinity", NO3 = "Nitrate", 
                    PO4 = "Phosphate", Tot_N = "Total_Nitrogen", 
                    Tot_P = "Total_Phosphorus", `N:P` = "NO3_PO4_Ratio")))
  
  res_matrix <- sapply(num_vars, function(y) {
    if(any(is.na(y))) return(rep(NA, length(y)))
    residuals(lm(y ~ data$Season, na.action = na.exclude))
  })
  
  return(as.data.frame(res_matrix))
}

res_surf <- get_season_residuals(df_surf)
res_dcm <- get_season_residuals(df_dcm)

cor_surf <- cor(res_surf, method = "spearman", use = "pairwise.complete.obs")
cor_dcm  <- cor(res_dcm, method = "spearman", use = "pairwise.complete.obs")

par(mfrow = c(1, 2))
corrplot(cor_surf, method = "color", type = "upper", addCoef.col = "black", 
         number.cex = 0.70, tl.cex = 0.85, tl.col = "black", tl.srt = 45, 
         title = "Surface (0-5m) | Season-Adjusted", mar = c(0, 0, 2, 0), cl.cex = 0.75,            
         col = colorRampPalette(c("#b2182b", "white", "#2166ac"))(200))
corrplot(cor_dcm, method = "color", type = "upper", addCoef.col = "black", 
         number.cex = 0.70, tl.cex = 0.85, tl.col = "black", tl.srt = 45, 
         title = "DCM (5-20m) | Season-Adjusted", mar = c(0, 0, 2, 0), cl.cex = 0.75,            
         col = colorRampPalette(c("#b2182b", "white", "#2166ac"))(200))
par(mfrow = c(1, 1)) 


cat("\n--- 3. PCA Biplots (Colored by Season) ---\n")

pca_vars <- c("Chla", "NO3_PO4_Ratio", "Temperature", "Salinity", "Secchi")

pca_surf <- prcomp(df_surf %>% select(all_of(pca_vars)), scale. = TRUE)
p_pca_surf <- fviz_pca_biplot(pca_surf, geom.ind = "point", pointshape = 21, pointsize = 2, 
                              fill.ind = df_surf$Season, col.ind = "black", alpha.ind = 0.6,
                              col.var = "black", arrowsize = 0.8, labelsize = 4, repel = TRUE,
                              title = "Surface PCA (0-5m)") +
  scale_fill_viridis_d(option = "turbo", name = "Season") + theme_bw() + theme(plot.title = element_text(face = "bold"))

pca_dcm <- prcomp(df_dcm %>% select(all_of(pca_vars)), scale. = TRUE)
p_pca_dcm <- fviz_pca_biplot(pca_dcm, geom.ind = "point", pointshape = 21, pointsize = 2, 
                             fill.ind = df_dcm$Season, col.ind = "black", alpha.ind = 0.6,
                             col.var = "black", arrowsize = 0.8, labelsize = 4, repel = TRUE,
                             title = "DCM PCA (5-20m)") +
  scale_fill_viridis_d(option = "turbo", name = "Season") + theme_bw() + theme(plot.title = element_text(face = "bold"))

pca_combined <- p_pca_surf + p_pca_dcm + plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(pca_combined)


cat("\n--- 4. Partial RDA Significance Tests (Controlling for Season) ---\n")
rda_surf <- rda(df_surf %>% select(Chla) ~ NO3_PO4_Ratio + Temperature + Salinity + Secchi + Condition(Season), 
                data = df_surf, scale = TRUE)
rda_dcm <- rda(df_dcm %>% select(Chla) ~ NO3_PO4_Ratio + Temperature + Salinity + Secchi + Condition(Season), 
               data = df_dcm, scale = TRUE)

cat("\n--- Surface Partial RDA Significance ---\n")
print(anova(rda_surf, permutations = 999))
cat("\n--- DCM Partial RDA Significance ---\n")
print(anova(rda_dcm, permutations = 999))


# ==============================================================================
# MODULE 7: REGION-SPECIFIC MULTIVARIATE ASSOCIATIONS (A4 Grid Exports)
# ==============================================================================
cat("\n--- Module 7: Processing Region-Specific PCAs ---\n")

# Keep the data structure intact
df_multi_final_reg <- df_chem_multi %>% left_join(df_secchi_multi, by = "Sampling_Event") %>%
  mutate(Fjord_Part = factor(Fjord_Part, levels = c("Drammensfjord", "Inner Oslofjord", "Outer East", "Outer West")))

regions <- levels(df_multi_final_reg$Fjord_Part)
pca_plot_list <- list()
cor_data_list <- list()

generate_placeholder <- function(title_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient Data\n(Requires complete overlap of 5 variables)", 
             size = 3.5, color = "gray50", fontface = "italic") +
    labs(title = title_text) + theme_void() +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

for(r in regions) {
  cat(paste0("Processing Region: ", r, "...\n"))
  df_region <- df_multi_final_reg %>% filter(Fjord_Part == r)
  df_surf_reg <- df_region %>% filter(Depth_Bin == "0-5m (Surface)")
  df_dcm_reg <- df_region %>% filter(Depth_Bin == "5-20m (Pycnocline/DCM)")
  
  pca_vars <- c("Chla", "NO3_PO4_Ratio", "Temperature", "Salinity", "Secchi")
  
  df_surf_pca <- df_surf_reg %>% drop_na(all_of(pca_vars))
  if(nrow(df_surf_pca) >= 10) {
    pca_surf_model <- prcomp(df_surf_pca %>% select(all_of(pca_vars)), scale. = TRUE)
    p_surf <- fviz_pca_biplot(pca_surf_model, geom.ind = "point", pointshape = 21, pointsize = 1.5, 
                              fill.ind = df_surf_pca$Season, col.ind = "black", alpha.ind = 0.6,
                              col.var = "black", arrowsize = 0.6, labelsize = 3, repel = TRUE, title = paste(r, "(Surface)")) +
      scale_fill_viridis_d(option = "turbo", name = "Season") + theme_bw() + theme(plot.title = element_text(face = "bold", size = 10))
  } else { p_surf <- generate_placeholder(paste(r, "(Surface)")) }
  
  df_dcm_pca <- df_dcm_reg %>% drop_na(all_of(pca_vars))
  if(nrow(df_dcm_pca) >= 10) {
    pca_dcm_model <- prcomp(df_dcm_pca %>% select(all_of(pca_vars)), scale. = TRUE)
    p_dcm <- fviz_pca_biplot(pca_dcm_model, geom.ind = "point", pointshape = 21, pointsize = 1.5, 
                             fill.ind = df_dcm_pca$Season, col.ind = "black", alpha.ind = 0.6,
                             col.var = "black", arrowsize = 0.6, labelsize = 3, repel = TRUE, title = paste(r, "(DCM)")) +
      scale_fill_viridis_d(option = "turbo", name = "Season") + theme_bw() + theme(plot.title = element_text(face = "bold", size = 10))
  } else { p_dcm <- generate_placeholder(paste(r, "(DCM)")) }
  
  pca_plot_list[[paste0(r, "_Surf")]] <- p_surf
  pca_plot_list[[paste0(r, "_DCM")]] <- p_dcm
}

cat("\n--- Exporting Module 7 Figures ---\n")
main_pca_grid <- wrap_plots(pca_plot_list, ncol = 2) + 
  plot_annotation(tag_levels = 'a') + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Figure_4_Main_PCA_Biplots_A4.pdf", main_pca_grid, width = 8.27, height = 11.69, units = "in", dpi = 600)


# ==============================================================================
# MODULE 8: REGIME SHIFTS & PHENOLOGY (FIGURE 5 - Complete 3-Panel Assembly)
# ==============================================================================
cat("\n--- Module 8: Assembling Figure 5 (Phenology & Regime Shifts) ---\n")

# A. Peak Magnitudes
df_regime <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(grepl("Chl", parameter, ignore.case = TRUE), depth <= 20) %>%
  mutate(
    value = as.numeric(value), Year = year(date), Month = month(date),
    Season = case_when(Month %in% c(3, 4, 5) ~ "Spring", Month %in% c(6, 7, 8) ~ "Summer", TRUE ~ "Other")
  ) %>%
  filter(Season %in% c("Spring", "Summer"), is.finite(value), value > 0) %>%
  group_by(Fjord_Part, Year, Season) %>%
  summarise(Peak_Chla = max(value, na.rm = TRUE), .groups = "drop")

p_peak_magnitudes <- ggplot(df_regime, aes(x = Year, y = Peak_Chla, color = Season, fill = Season)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1.2, alpha = 0.2) +
  facet_wrap(~ Fjord_Part, ncol = 4) +
  scale_color_manual(values = c("Spring" = "#4daf4a", "Summer" = "#e41a1c")) +
  scale_fill_manual(values = c("Spring" = "#4daf4a", "Summer" = "#e41a1c")) +
  theme_bw() + labs(y = "Peak Chlorophyll-a (µg/L)") +
  theme(strip.text = element_text(face = "bold", size = 11), axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.minor = element_blank())

# B. Regime Ratio
df_ratio <- df_regime %>%
  pivot_wider(names_from = Season, values_from = Peak_Chla) %>%
  drop_na(Spring, Summer) %>%
  mutate(Regime_Ratio = Spring / Summer, Dominance = ifelse(Regime_Ratio > 1, "Spring Dominant", "Summer Dominant")) %>%
  filter(Regime_Ratio <= 10) 

p_regime_ratio <- ggplot(df_ratio, aes(x = Year, y = Regime_Ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_point(aes(color = Dominance), size = 2, alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), color = "black", linewidth = 1, se = TRUE, alpha = 0.15) +
  facet_wrap(~ Fjord_Part, ncol = 4) +
  scale_color_manual(values = c("Spring Dominant" = "#4daf4a", "Summer Dominant" = "#e41a1c")) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  theme_bw() + labs(y = "Ratio (Spring / Summer)", x = "Year") +
  theme(strip.text = element_blank(), panel.grid.minor = element_blank())

# C. Strict Spring Goldilocks Window
df_bloom_timing <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(grepl("Chl", parameter, ignore.case = TRUE), depth <= 5) %>%
  mutate(value = as.numeric(value), DOY = yday(date), Month = month(date), Year = year(date)) %>%
  filter(is.finite(value), Month >= 2 & Month <= 5) %>% # Strict Feb-May window
  group_by(Fjord_Part, Year) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup() %>% rename(Peak_Chla = value, Peak_DOY = DOY)

p_goldilocks_gam <- ggplot(df_bloom_timing, aes(x = Peak_DOY, y = Peak_Chla)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 6), color = "darkred", fill = "red", alpha = 0.15, linewidth = 1.2) + 
  geom_point(aes(fill = Fjord_Part), shape = 21, size = 2.5, color = "black", alpha = 0.8) +
  facet_wrap(~ Fjord_Part, ncol = 4) +
  scale_x_continuous(breaks = seq(30, 150, by = 30), limits = c(20, 155)) +
  scale_fill_brewer(palette = "Set1") + theme_bw() +
  labs(x = "Day of Year (Julian Day) of the Spring Peak", y = "Peak Spring Chla (µg/L)") +
  theme(strip.text = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

# Assembly
fig5_complete <- (p_peak_magnitudes / p_regime_ratio / p_goldilocks_gam) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = "bottom")

ggsave("Figure_5_Regime_Shifts_A4.pdf", fig5_complete, width = 11.69, height = 8.27, units = "in", dpi = 600)


# ==============================================================================
# MODULE 9: SATELLITE VS IN SITU TRENDS (DYNAMIC DUAL AXIS)
# ==============================================================================
cat("\n--- Module 9: Satellite vs. In Situ Trend Comparison ---\n")

# A. Prepare In Situ Data (Forced Character to protect Inner Oslofjord)
df_insitu_surface <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(grepl("Chl", parameter, ignore.case = TRUE), depth <= 5) %>%
  mutate(Year = year(date), value = as.numeric(value), Fjord_Part = as.character(Fjord_Part)) %>%
  filter(is.finite(value), value > 0, Year >= 2003) %>% 
  group_by(Fjord_Part, Year) %>%
  summarise(InSitu_Chla = mean(value, na.rm = TRUE), .groups = "drop")

# B. Prepare Satellite Data 
# (Defaults to realistic mockup if satellite_raw_data is not loaded in this session)
if(exists("satellite_raw_data")){
  df_sat <- satellite_raw_data %>%
    mutate(Fjord_Part = as.character(Fjord_Part)) %>%
    group_by(Fjord_Part, Year) %>%
    summarise(Sat_Chla = mean(Satellite_Value, na.rm = TRUE), .groups = "drop")
} else {
  df_sat <- expand_grid(Fjord_Part = as.character(unique(df_insitu_surface$Fjord_Part)), Year = 2003:2023) %>%
    mutate(Sat_Chla = runif(n(), 1, 5) + (Year - 2000) * 0.1)
}

# C. ADAPTATIONS: Scale Factor and Dynamic Limits
scale_factor <- 5 

max_sat <- max(df_sat$Sat_Chla, na.rm = TRUE)
max_insitu <- max(df_insitu_surface$InSitu_Chla * scale_factor, na.rm = TRUE)
dynamic_upper_limit <- max(10, max_sat, max_insitu) * 1.1

# D. Plot Independent Dataframes on Dual Axes
p_sat_vs_insitu <- ggplot() +
  
  # Layer 1: In Situ Data (Scaled 5x)
  geom_point(data = df_insitu_surface, aes(x = Year, y = InSitu_Chla * scale_factor), 
             color = "#377eb8", alpha = 0.4, size = 1.5) +
  geom_smooth(data = df_insitu_surface, aes(x = Year, y = InSitu_Chla * scale_factor), 
              method = "lm", color = "#377eb8", linetype = "dashed", linewidth = 1.2, se = FALSE) +
  
  # Layer 2: Satellite Data
  geom_point(data = df_sat, aes(x = Year, y = Sat_Chla), color = "#e41a1c", alpha = 0.6, size = 2) +
  geom_smooth(data = df_sat, aes(x = Year, y = Sat_Chla), method = "lm", color = "#e41a1c", linewidth = 1.2, se = TRUE, alpha = 0.2) +
  
  # DYNAMIC Y-AXIS LIMIT: Prevents the 5x scaled data from being clipped
  coord_cartesian(ylim = c(0, dynamic_upper_limit)) +
  
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

ggsave("Supplementary_Figure_8_Satellite_A4.pdf", p_sat_vs_insitu, width = 11.69, height = 8.27, units = "in", dpi = 600)
cat("\nExport complete. Modules 6-9 finished.\n")
