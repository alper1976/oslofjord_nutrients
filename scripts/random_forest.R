# ==============================================================================
# MODULE 11: CONDITIONAL RANDOM FOREST (Starting from Master Object)
# ==============================================================================
library(party)
library(caret)
library(Metrics)
library(broom)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(ggplot2)

cat("\n--- 1. Preparing Data directly from oslo_final_sf ---\n")

# A. Extract Master Data and Force Temporal Columns
df_master <- oslo_final_sf %>%
  st_drop_geometry() %>%
  mutate(
    date = as.Date(date),
    Year = year(date),
    Month = month(date),
    # Create the ultimate, foolproof grouping ID for cross-validation!
    # All depths sampled in the same Fjord Part on the exact same day get the same ID.
    Sampling_Event = as.factor(paste(Fjord_Part, date, sep = "_"))
  ) %>%
  filter(!is.na(date))

# B. Extract Secchi Depth
df_secchi <- df_master %>%
  filter(grepl("Secchi", parameter, ignore.case = TRUE)) %>%
  mutate(value = as.numeric(value)) %>%
  filter(is.finite(value), value > 0) %>%
  group_by(Sampling_Event) %>% # Match directly by the unique event
  summarise(Secchi = mean(value, na.rm = TRUE), .groups = "drop")

# C. Extract Chemistry & Chlorophyll
df_chem <- df_master %>%
  filter(!grepl("Secchi", parameter, ignore.case = TRUE)) %>%
  mutate(
    parameter_clean = case_when(
      grepl("Chl", parameter, ignore.case = TRUE) ~ "Chla",
      grepl("Nitrate", parameter, ignore.case = TRUE) ~ "Nitrate",
      grepl("Phosphate", parameter, ignore.case = TRUE) ~ "Phosphate",
      grepl("Temp", parameter, ignore.case = TRUE) ~ "Temperature",
      grepl("Salin", parameter, ignore.case = TRUE) ~ "Salinity",
      TRUE ~ "IGNORE"
    ),
    Depth_Bin = case_when(
      depth <= 5 ~ "0-5m (Surface)",
      depth > 5 & depth <= 20 ~ "5-20m (Pycnocline/DCM)",
      TRUE ~ "Deep"
    ),
    # Convert nutrients to molarity for accurate stoichiometry
    value_clean = case_when(
      parameter_clean == "Nitrate" ~ as.numeric(value) / 14.01,
      parameter_clean == "Phosphate" ~ as.numeric(value) / 30.97,
      TRUE ~ as.numeric(value)
    )
  ) %>%
  filter(parameter_clean != "IGNORE", Depth_Bin != "Deep", is.finite(value_clean)) %>%
  group_by(Sampling_Event, Fjord_Part, Depth_Bin, Year, Month, parameter_clean) %>%
  summarise(value = mean(value_clean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter_clean, values_from = value) %>%
  mutate(NO3_PO4_Ratio = Nitrate / (Phosphate + 0.01)) # Avoid dividing by zero

# D. Merge and Finalize Model Data
rf_data <- df_chem %>%
  inner_join(df_secchi, by = "Sampling_Event") %>%
  drop_na(Chla, Temperature, Salinity, Secchi, NO3_PO4_Ratio) %>%
  mutate(
    Fjord_Part = as.factor(Fjord_Part),
    Depth_Bin = as.factor(Depth_Bin)
  ) %>%
  select(Chla, Temperature, Salinity, Secchi, NO3_PO4_Ratio, 
         Fjord_Part, Depth_Bin, Sampling_Event) %>%
  rename(Target = Chla)

cat("\n--- 2. Iterative Random Forest Setup ---\n")

set.seed(1000)
max_num_of_iterations <- 100

res <- NULL; pred <- NULL; true <- NULL
rmse_list <- NULL; mae_list <- NULL; imp <- NULL
rff <- list()

pb <- txtProgressBar(min = 0, max = max_num_of_iterations, style = 3, width = 100, char = "=")

# 3. The Loop
for(i in 1:max_num_of_iterations){  
  
  # Grouped Split: Keep exact dates/locations together!
  events <- unique(rf_data$Sampling_Event)
  train_events <- sample(events, size = 0.7 * length(events))
  
  training <- rf_data %>% filter(Sampling_Event %in% train_events)
  testing  <- rf_data %>% filter(!(Sampling_Event %in% train_events))
  
  train_model_data <- training %>% select(-Sampling_Event)
  test_model_data  <- testing %>% select(-Sampling_Event)
  
  folds <- groupKFold(training$Sampling_Event, k = 10)
  repeat_cv <- trainControl(method = 'cv', index = folds, savePredictions = "final")
  
  suppressWarnings( 
    model <- train(Target ~ ., data = train_model_data, 
                   method = 'cforest', 
                   trControl = repeat_cv,
                   controls = cforest_unbiased(ntree = 500), 
                   tuneGrid = expand.grid(mtry = 3)) 
  )
  
  rff[[i]] <- model
  p <- predict(model, test_model_data)
  t <- test_model_data$Target
  
  imp1 <- t(as.data.frame(varimp(model$finalModel)))
  imp <- rbind(imp, imp1)
  
  p.lm <- lm(p ~ t)
  p.lm.slope <- summary(p.lm)
  res1 <- glance(p.lm) %>% mutate(slo = p.lm.slope$coefficients[2], ints = p.lm.slope$coefficients[1])
  res <- rbind(res, res1) 
  
  pred <- rbind(pred, data.frame(p = p, Iteration = i)) 
  true <- rbind(true, data.frame(t = t, Iteration = i)) 
  rmse_list <- rbind(rmse_list, rmse(t, p))
  
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n--- 3. Random Forest Training Complete ---\n")

# 4. Plot Conditional Variable Importance
varim_mean <- colMeans(imp)
varim_df <- data.frame(Variable = names(varim_mean), MSE = varim_mean) %>%
  arrange(MSE) %>%
  mutate(Variable = factor(Variable, levels = Variable))

p_rf_importance <- ggplot(data = varim_df, aes(x = Variable, y = MSE)) +
  geom_point(stat = 'identity', size = 4, color = "#2166ac") +
  geom_linerange(aes(ymin = 0, ymax = MSE), color = "#2166ac", linewidth = 1) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Primary Drivers of Phytoplankton Biomass (Chla)",
    subtitle = "Conditional Random Forest Variable Importance",
    x = "Predictors", y = "Importance (Increase in MSE when removed)"
  )

print(p_rf_importance)
