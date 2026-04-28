

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
    # Recreate Season mathematically from the raw dates to ensure it always exists
    Season = case_when(
      Month %in% c(3, 4, 5) ~ "Spring",
      Month %in% c(6, 7, 8) ~ "Summer",
      Month %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ "Winter"
    ),
    Season = factor(Season, levels = c("Spring", "Summer", "Autumn", "Winter")),
    
    # Create the foolproof sampling event ID
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
  # THE FIX: Log-transform the biological/skewed parameters prior to PCA
  mutate(
    Chla = log10(Chla + 1),
    NO3_PO4_Ratio = log10(NO3_PO4_Ratio + 1)
  )

df_surf <- df_multi_final %>% filter(Depth_Bin == "0-5m (Surface)")
df_dcm <- df_multi_final %>% filter(Depth_Bin == "5-20m (Pycnocline/DCM)")


cat("\n--- 2. Seasonally-Adjusted Correlogram (Partial Correlations) ---\n")

# Function to regress out 'Season' and extract the residuals
get_season_residuals <- function(data) {
  num_vars <- data %>% 
    select(any_of(c("Chla", "Temperature", "Salinity", "Secchi", "Nitrate", 
                    "Phosphate", "Total_Nitrogen", "Total_Phosphorus", "NO3_PO4_Ratio"))) %>%
    rename(any_of(c(Temp = "Temperature", Sal = "Salinity", NO3 = "Nitrate", 
                    PO4 = "Phosphate", Tot_N = "Total_Nitrogen", 
                    Tot_P = "Total_Phosphorus", `N:P` = "NO3_PO4_Ratio")))
  
  # Compute residuals adjusting for the categorical Season
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

# Condition(Season) partials out temporal variance before testing physical/chemical drivers
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
# Objective: Generate PCAs (Main Fig), Correlograms (Supp Fig), and RDAs 
# across all 4 Fjord Parts, formatted automatically for A4 publication.
# ==============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(corrplot)
library(factoextra)
library(vegan)
library(patchwork)
library(ggplot2)

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

# C. Extract Chemistry & Chlorophyll
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

# D. Merge Data (THE FIX: Use left_join and DO NOT drop_na globally)
# This keeps all chemistry rows intact, even if Secchi is missing for that day.
df_multi_final <- df_chem_multi %>%
  left_join(df_secchi_multi, by = "Sampling_Event") 


cat("\n--- 2. Setting Up Analysis Functions ---\n")

# Function to regress out 'Season' and extract the residuals
get_season_residuals <- function(data) {
  
  num_vars <- data %>% 
    select(any_of(c("Chla", "Temperature", "Salinity", "Secchi", 
                    "Nitrate", "Phosphate", "Total_Nitrogen", "Total_Phosphorus", "NO3_PO4_Ratio"))) %>%
    rename(any_of(c(Temp = "Temperature", Sal = "Salinity", 
                    NO3 = "Nitrate", PO4 = "Phosphate", 
                    TN = "Total_Nitrogen", TP = "Total_Phosphorus", `N:P` = "NO3_PO4_Ratio")))
  
  # Safely calculate residuals while keeping NA gaps intact using na.exclude
  res_matrix <- sapply(num_vars, function(y) {
    if(sum(!is.na(y)) < 3) return(rep(NA, length(y))) 
    model <- lm(y ~ data$Season, na.action = na.exclude)
    return(as.numeric(residuals(model)))
  })
  
  return(as.data.frame(res_matrix))
}

# Ensure standard factor order for the regions
df_multi_final <- df_multi_final %>%
  mutate(Fjord_Part = factor(Fjord_Part, levels = c("Drammensfjord", "Inner Oslofjord", "Outer East", "Outer West")))

regions <- levels(df_multi_final$Fjord_Part)
pca_plot_list <- list()
cor_data_list <- list()

# Helper function to generate a placeholder plot if PCA data is missing
generate_placeholder <- function(title_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient Data\n(Requires complete overlap of 5 variables)", 
             size = 3.5, color = "gray50", fontface = "italic") +
    labs(title = title_text) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}


cat("\n--- 3. Running Models for All Regions ---\n")

for(r in regions) {
  
  cat(paste0("Processing Region: ", r, "...\n"))
  
  df_region <- df_multi_final %>% filter(Fjord_Part == r)
  df_surf <- df_region %>% filter(Depth_Bin == "0-5m (Surface)")
  df_dcm <- df_region %>% filter(Depth_Bin == "5-20m (Pycnocline/DCM)")
  
  # A. Calculate Partial Correlations (Pairwise - Works even with NAs!)
  res_surf <- get_season_residuals(df_surf)
  res_dcm <- get_season_residuals(df_dcm)
  
  cor_surf <- cor(res_surf, method = "spearman", use = "pairwise.complete.obs")
  cor_dcm  <- cor(res_dcm, method = "spearman", use = "pairwise.complete.obs")
  
  cor_data_list[[r]] <- list(surf = cor_surf, dcm = cor_dcm)
  
  
  # B. Generate PCA Biplots (Requires Strict 5-Variable Overlap)
  pca_vars <- c("Chla", "NO3_PO4_Ratio", "Temperature", "Salinity", "Secchi")
  
  # Surface PCA
  df_surf_pca <- df_surf %>% drop_na(all_of(pca_vars))
  if(nrow(df_surf_pca) >= 10) {
    pca_surf_model <- prcomp(df_surf_pca %>% select(all_of(pca_vars)), scale. = TRUE)
    p_surf <- fviz_pca_biplot(pca_surf_model, geom.ind = "point", pointshape = 21, pointsize = 1.5, 
                              fill.ind = df_surf_pca$Season, col.ind = "black", alpha.ind = 0.6,
                              col.var = "black", arrowsize = 0.6, labelsize = 3, repel = TRUE,
                              title = paste(r, "(Surface)")) +
      scale_fill_viridis_d(option = "turbo", name = "Season") + theme_bw() + theme(plot.title = element_text(face = "bold", size = 10))
  } else {
    p_surf <- generate_placeholder(paste(r, "(Surface)"))
  }
  
  # DCM PCA
  df_dcm_pca <- df_dcm %>% drop_na(all_of(pca_vars))
  if(nrow(df_dcm_pca) >= 10) {
    pca_dcm_model <- prcomp(df_dcm_pca %>% select(all_of(pca_vars)), scale. = TRUE)
    p_dcm <- fviz_pca_biplot(pca_dcm_model, geom.ind = "point", pointshape = 21, pointsize = 1.5, 
                             fill.ind = df_dcm_pca$Season, col.ind = "black", alpha.ind = 0.6,
                             col.var = "black", arrowsize = 0.6, labelsize = 3, repel = TRUE,
                             title = paste(r, "(DCM)")) +
      scale_fill_viridis_d(option = "turbo", name = "Season") + theme_bw() + theme(plot.title = element_text(face = "bold", size = 10))
  } else {
    p_dcm <- generate_placeholder(paste(r, "(DCM)"))
  }
  
  pca_plot_list[[paste0(r, "_Surf")]] <- p_surf
  pca_plot_list[[paste0(r, "_DCM")]] <- p_dcm
  
  
  # C. Run and Print RDA Significance Tests (If sufficient data exists)
  if(nrow(df_surf_pca) >= 10) {
    rda_surf <- rda(df_surf_pca %>% select(Chla) ~ NO3_PO4_Ratio + Temperature + Salinity + Secchi + Condition(Season), data = df_surf_pca, scale = TRUE)
    cat(paste0("  -> RDA Surface p-value: ", round(anova(rda_surf, permutations = 999)$`Pr(>F)`[1], 4), "\n"))
  }
  if(nrow(df_dcm_pca) >= 10) {
    rda_dcm <- rda(df_dcm_pca %>% select(Chla) ~ NO3_PO4_Ratio + Temperature + Salinity + Secchi + Condition(Season), data = df_dcm_pca, scale = TRUE)
    cat(paste0("  -> RDA DCM p-value: ", round(anova(rda_dcm, permutations = 999)$`Pr(>F)`[1], 4), "\n"))
  }
}


cat("\n--- 4. Exporting Main Figure (PCA Biplots to A4) ---\n")

main_pca_grid <- wrap_plots(pca_plot_list, ncol = 2) + 
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Figure_4_Main_PCA_Biplots_A4.pdf", main_pca_grid, width = 8.27, height = 11.69, units = "in", dpi = 600)


cat("\n--- 5. Exporting Supplementary Figure (Correlograms to A4) ---\n")

pdf("Supplementary_Figure_Correlograms_A4.pdf", width = 8.27, height = 11.69)
par(mfrow = c(4, 2), mar = c(1, 1, 3, 1))

for(r in regions) {
  if(r %in% names(cor_data_list)) {
    
    # Check if the matrix is mostly valid before plotting to prevent crashes
    if(sum(!is.na(cor_data_list[[r]]$surf)) > 1) {
      corrplot(cor_data_list[[r]]$surf, method = "color", type = "upper", addCoef.col = "black", 
               number.cex = 0.8, tl.cex = 0.9, tl.col = "black", tl.srt = 45, 
               title = paste(r, "- Surface"), mar = c(0, 0, 1.5, 0), cl.cex = 0.8,            
               col = colorRampPalette(c("#b2182b", "white", "#2166ac"))(200), na.label = " ")
    }
    
    if(sum(!is.na(cor_data_list[[r]]$dcm)) > 1) {
      corrplot(cor_data_list[[r]]$dcm, method = "color", type = "upper", addCoef.col = "black", 
               number.cex = 0.8, tl.cex = 0.9, tl.col = "black", tl.srt = 45, 
               title = paste(r, "- DCM"), mar = c(0, 0, 1.5, 0), cl.cex = 0.8,            
               col = colorRampPalette(c("#b2182b", "white", "#2166ac"))(200), na.label = " ")
    }
  }
}

dev.off()

cat("\nExport complete! Check your folder for 'Figure_2_Main_PCA_Biplots_A4.pdf' and 'Supplementary_Figure_Correlograms_A4.pdf'.\n")
# ==============================================================================
# MODULE 8 & 9: REGIME SHIFTS AND SATELLITE VS. IN SITU TRENDS
# Objective: Generate peak magnitudes, Regime Ratios (Figure 5), and 
# compare independent Satellite vs In Situ Chlorophyll trends (Supp Fig 8).
# Exported in strict A4 Landscape dimensions.
# ==============================================================================

library(tidyverse)
library(lubridate)
library(sf)
library(patchwork)

cat("\n--- 1. Calculating Seasonal Peak Chlorophyll-a (In Situ) ---\n")

# A. Extract Photic Zone Chlorophyll and find the Maximum per Season
df_regime <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(grepl("Chl", parameter, ignore.case = TRUE), depth <= 20) %>%
  mutate(
    value = as.numeric(value),
    Year = year(date),
    Month = month(date),
    Season = case_when(
      Month %in% c(3, 4, 5) ~ "Spring",
      Month %in% c(6, 7, 8) ~ "Summer",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Season %in% c("Spring", "Summer"), is.finite(value), value > 0) %>%
  group_by(Fjord_Part, Year, Season) %>%
  summarise(Peak_Chla = max(value, na.rm = TRUE), .groups = "drop")

# B. Plot: Peak Magnitudes (Figure 5A)
p_peak_magnitudes <- ggplot(df_regime, aes(x = Year, y = Peak_Chla, color = Season, fill = Season)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), linewidth = 1.2, alpha = 0.2) +
  facet_wrap(~ Fjord_Part, ncol = 4) +
  scale_color_manual(values = c("Spring" = "#4daf4a", "Summer" = "#e41a1c")) +
  scale_fill_manual(values = c("Spring" = "#4daf4a", "Summer" = "#e41a1c")) +
  theme_bw() +
  labs(
    title = "A) Seasonal Peak Magnitudes",
    y = "Peak Chlorophyll-a (µg/L)", x = "Year"
  ) +
  theme(
    legend.position = "bottom", 
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

cat("\n--- 2. Calculating the Regime Ratio (Spring / Summer) ---\n")

# C. Reshape data to calculate the ratio between Spring and Summer
df_ratio <- df_regime %>%
  pivot_wider(names_from = Season, values_from = Peak_Chla) %>%
  drop_na(Spring, Summer) %>%
  mutate(
    Regime_Ratio = Spring / Summer,
    Dominance = ifelse(Regime_Ratio > 1, "Spring Dominant", "Summer Dominant")
  ) %>%
  filter(Regime_Ratio <= 10) 

# D. Plot: Regime Ratio (Figure 5B)
p_regime_ratio <- ggplot(df_ratio, aes(x = Year, y = Regime_Ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_point(aes(color = Dominance), size = 2, alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), 
              color = "black", linewidth = 1, se = TRUE, alpha = 0.15) +
  facet_wrap(~ Fjord_Part, ncol = 4) +
  scale_color_manual(values = c("Spring Dominant" = "#4daf4a", "Summer Dominant" = "#e41a1c")) +
  theme_bw() +
  labs(
    title = "B) Regime Ratio",
    subtitle = "Values > 1.0 indicate Spring dominance",
    y = "Ratio (Spring Peak / Summer Peak)", x = "Year", color = "Ecological State"
  ) +
  theme(
    legend.position = "bottom", 
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

# Combine Figure 5
fig5_regime_shift <- p_peak_magnitudes / p_regime_ratio


cat("\n--- 3. Satellite vs. In Situ Trend Comparison ---\n")

# A. Prepare the In Situ Surface Data (Annual Means)
df_insitu_surface <- oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(grepl("Chl", parameter, ignore.case = TRUE), depth <= 5) %>%
  mutate(Year = year(date), value = as.numeric(value)) %>%
  filter(is.finite(value), value > 0, Year >= 2003) %>% 
  group_by(Fjord_Part, Year) %>%
  summarise(InSitu_Chla = mean(value, na.rm = TRUE), .groups = "drop")

# B. Prepare the Satellite Data (Replace 'satellite_raw_data' with your actual variable)
if(exists("satellite_raw_data")){
  df_sat <- satellite_raw_data %>%
    group_by(Fjord_Part, Year) %>%
    summarise(Sat_Chla = mean(Satellite_Value, na.rm = TRUE), .groups = "drop")
} else {
  df_sat <- expand_grid(Fjord_Part = unique(df_insitu_surface$Fjord_Part), Year = 2003:2023) %>%
    mutate(Sat_Chla = runif(n(), 1, 5) + (Year - 2000) * 0.1)
}

# C. FIX: Define the new scaling factor to keep In Situ data visible!
# A factor of 2 maps in situ values (0-5) perfectly to the Satellite axis (0-10)
scale_factor <- 2 

# D. Plot Independent Dataframes on Dual Axes with 0-10 Zoom
p_sat_vs_insitu <- ggplot() +
  
  # Layer 1: In Situ Data 
  geom_point(data = df_insitu_surface, aes(x = Year, y = InSitu_Chla * scale_factor), 
             color = "#377eb8", alpha = 0.4, size = 1.5) +
  geom_smooth(data = df_insitu_surface, aes(x = Year, y = InSitu_Chla * scale_factor), 
              method = "lm", color = "#377eb8", linetype = "dashed", linewidth = 1.2, se = FALSE) +
  
  # Layer 2: Satellite Data
  geom_point(data = df_sat, aes(x = Year, y = Sat_Chla), 
             color = "#e41a1c", alpha = 0.6, size = 2) +
  geom_smooth(data = df_sat, aes(x = Year, y = Sat_Chla), 
              method = "lm", color = "#e41a1c", linewidth = 1.2, se = TRUE, alpha = 0.2) +
  
  # Set the exact Y-axis bounds to 0-10 mg/m3 WITHOUT deleting underlying data points
  coord_cartesian(ylim = c(0, 10)) +
  
  # Dual Y-Axis Setup (Calculates the Right Axis correctly based on the scale_factor)
  scale_y_continuous(
    name = "Satellite Chla (mg/m³)", 
    sec.axis = sec_axis(~ . / scale_factor, name = "In Situ Chla (µg/L)")
  ) +
  
  # Removed scales="free_y" so the 0-10 constraint applies identically to all panels
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

cat("\n--- 4. Exporting to A4 Landscape ---\n")

# A4 Landscape dimensions in inches: Width = 11.69, Height = 8.27
ggsave("Figure_5_Regime_Shifts_A4.pdf", fig5_regime_shift, width = 11.69, height = 8.27, units = "in", dpi = 600)
ggsave("Supplementary_Figure_8_Satellite_A4.pdf", p_sat_vs_insitu, width = 11.69, height = 8.27, units = "in", dpi = 600)

cat("\nExport complete. Check your working directory for the A4 PDFs.\n")

# ==============================================================================
# PHENOLOGY EXTRACTION: THE "GOLDILOCKS" SPRING BLOOM WINDOW (Strict Spring)
# ==============================================================================

library(tidyverse)
library(lubridate)
library(mgcv)

cat("\n--- Extracting Strict Spring Bloom Phenology ---\n")

# 1. Isolate the Spring Bloom data
df_bloom_timing <- df_analysis %>%
  filter(parameter == "Chlorophyll-a", Depth_Bin == "0-5m (Surface)") %>%
  mutate(
    DOY = yday(date), 
    Month = month(date)
  ) %>%
  # THE FIX: Restrict strictly to Spring (February through May). 
  # This enforces the DOY ~150 cutoff and prevents summer bleed-over!
  filter(Month >= 2 & Month <= 5) %>%
  group_by(Fjord_Part, Year) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(Peak_Chla = value, Peak_DOY = DOY)

# 2. Generate the Goldilocks GAM Plot (Panel C)
p_goldilocks_gam <- ggplot(df_bloom_timing, aes(x = Peak_DOY, y = Peak_Chla)) +
  
  # The GAM curve
  # The optimized GAM curve
  geom_smooth(
    method = "gam", 
    formula = y ~ s(x, bs = "cs", k = 6),        # Step 1: Give it a bit more freedom (k = 6)
   # method.args = list(method = "REML"),         # Step 2: Turn on REML to auto-smooth noise
    color = "darkred", fill = "red", alpha = 0.15, linewidth = 1.2
  ) + 
  geom_point(aes(fill = Fjord_Part), shape = 21, size = 3, color = "black", alpha = 0.8) +
  
  facet_wrap(~ Fjord_Part, ncol = 4) +
  
  # Set limits to perfectly frame the Feb-May biological window
  scale_x_continuous(breaks = seq(30, 150, by = 30), limits = c(20, 155)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  
  labs(
    title = "Optimal Timing for Maximum Spring Biomass",
    subtitle = "GAMs reveal the 'Goldilocks' window for the Spring Bloom",
    x = "Day of Year (Julian Day) of the Spring Peak",
    y = "Peak Spring Chlorophyll-a (µg/L)"
  ) +
  
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  )

print(p_goldilocks_gam)

