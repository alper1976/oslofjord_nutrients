library(tidyverse)
library(sf)
library(leaflet)
library(dbscan)
library(ggspatial)
library(prettymapr)
library(ggspatial)
library(patchwork)
library(pals) # For the base qualitative colors


# --- 1. Prep for Clustering ---
# Transform to UTM (25833) so distances are in meters

# Prepare dataframe for Leaflet
leaflet_data <- oslo_final_sf %>%
  mutate(lng = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

leaflet(leaflet_data) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~lng, lat = ~lat,
    color = "black", weight = 1,
    fillColor = ~colorFactor("viridis", source)(source), # Using Viridis
    fillOpacity = 0.8,
    radius = 5,
    # Added 'source' and 'parameter' to the label for better insight
    label = ~paste0("Source: ", source, " | ", source, " (", parameter, ")")
  ) %>%
  addLegend(pal = colorFactor("viridis", leaflet_data$source), 
            values = ~source, title = "Source")
# 4. Extract Lat/Long back into columns if you need them for Leaflet
oslo_final_sf <- st_transform(oslo_final_sf, 4326)

# 2. Create the Map
# Calculate the range and create a sequence, then keep every other one
lon_bounds <- st_bbox(oslo_final_sf)
# Create a sequence of 0.1 or 0.2 degree intervals depending on the span
x_breaks <- seq(floor(lon_bounds["xmin"] * 5) / 5, 
                ceiling(lon_bounds["xmax"] * 5) / 5, 
                by = 0.1)
# Filter to keep every other tick
x_breaks_filtered <- x_breaks[seq(1, length(x_breaks), by = 2)]
ggplot(oslo_final_sf) +
  # Use 'osm' or 'cartolight' as they are currently the most reliable 
  # for the url_format required by ggspatial
  annotation_map_tile(type = "osm", zoom = 11, alpha = 0.8) + 
  
  # Data points with a thin black stroke for better contrast against map text
  layer_spatial(
    data = oslo_final_sf, 
    aes(color = source, fill = source), 
    size = 2.5, 
    shape = 21,    # Shape 21 allows for separate fill and stroke (border)
    stroke = 0.4,
    alpha = 0.8
  ) +
  
  # Scientific annotations
  annotation_scale(
    location = "tr", 
    width_hint = 0.4,
    pad_x = unit(0.38, "in"), 
    pad_y = unit(0.05, "in"), # Pushes it up
    text_cex = 0.7
  ) +
  annotation_north_arrow(
    location = "tr", 
    which_north = "true", 
    pad_x = unit(0.7, "in"), 
    pad_y = unit(0.2, "in"), # Pushes arrow even higher above the scale
    style = north_arrow_fancy_orienteering,
    height = unit(1, "cm"),
    width = unit(1, "cm")
  ) +
  
  # Publication-standard Colors
  scale_color_manual(values = c("Vannmiljø" = "black", "Archive" = "black")) + # Border colors
  scale_fill_brewer(palette = "Set1") + # Fill colors
  
  # X-axis: Every other tick
  scale_x_continuous(breaks = x_breaks_filtered) +
  
  # Theme and Legend positioning
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.background = element_rect(fill = alpha("white", 0.5)),
    axis.title = element_text(size = 9, color = "grey30"),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  
  labs(
    title = "Sampling locations: Oslofjord",
    fill = "Data Source",
    color = "Data Source", # Match color to fill for a clean legend
    x = "Longitude", 
    y = "Latitude",
    caption = "Base map: © OpenStreetMap contributors"
  ) +
  
  # Force vertical legend entries
  guides(fill = guide_legend(ncol = 1), color = "none")


# 1. Ensure site_cluster is a factor
oslo_final_sf$site_cluster <- as.factor(oslo_final_sf$site_cluster)
n_clust <- nlevels(oslo_final_sf$site_cluster)

# 2. Generate exactly 70 (or more) distinct colors
# We take a high-contrast palette and expand it to the required number
set.seed(42) 
base_colors <- as.vector(pals::polychrome()) # Get the base 36 distinct colors
custom_palette <- colorRampPalette(base_colors)(n_clust)
# Shuffle the palette so clusters 1, 2, 3 (neighbors) get very different colors
cluster_colors <- sample(custom_palette)

# 3. Define every other tick for X-axis
lon_range <- st_bbox(oslo_final_sf)
x_breaks <- seq(floor(lon_range["xmin"]*5)/5, ceiling(lon_range["xmax"]*5)/5, by = 0.1)
x_labels <- x_breaks[seq(1, length(x_breaks), by = 2)]

# 4. Plot
ggplot(oslo_final_sf) +
  annotation_map_tile(type = "cartolight", zoom = 11, alpha = 0.9) + 
  
  layer_spatial(
    data = oslo_final_sf, 
    aes(color = site_cluster, fill = site_cluster), 
    size = 2, 
    shape = 21, 
    stroke = 0.3
  ) +
  
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  
  # Set x-axis breaks (every other tick)
  scale_x_continuous(breaks = x_labels) +
  
  # Allow drawing outside the plot panel
  coord_sf(clip = "off") + 
  
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    # Create white space at the bottom (margin in points)
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 80, l = 10, unit = "pt"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15))
  ) +
  
  labs(
    title = "Oslofjord Sampling Clusters (1930-2025)",
    subtitle = paste(n_clust, "spatially distinct clusters (1km radius)"),
    x = "Longitude (°E)", 
    y = "Latitude (°N)"
  )
# Prepare the data: Count samples per Year and Cluster
cluster_summary <- oslo_final_sf %>%
  st_drop_geometry() %>% # Remove spatial overhead for faster processing
  count(Year, site_cluster) %>%
  # Ensure site_cluster is treated as a factor for the Y-axis
  mutate(site_cluster = factor(site_cluster, levels = rev(sort(unique(as.numeric(site_cluster))))))

# Create the Heatmap
ggplot(cluster_summary, aes(x = Year, y = site_cluster, fill = n)) +
  geom_tile(color = "white", linewidth = 0.1) +
  # Use a perceptually uniform color scale (viridis)
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Sample Count") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(cluster_summary$Year), 
                                  max(cluster_summary$Year), by = 5)) +
  labs(
    title = "Sampling Intensity by Cluster (1930-2025)",
    subtitle = "Number of observations per year for each 1km cluster",
    x = "Year",
    y = "Cluster ID"
  ) +
  theme(
    axis.text.y = element_text(size = 6), # Small text for 70 clusters
    panel.grid = element_blank(),
    legend.position = "bottom"
  )
# Find the top 10 clusters by total count
top_clusters <- oslo_final_sf %>%
  st_drop_geometry() %>%
  count(site_cluster, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  pull(site_cluster)

# Plot only those clusters
oslo_final_sf %>%
  st_drop_geometry() %>%
  filter(site_cluster %in% top_clusters) %>%
  ggplot(aes(x = Year, fill = source)) +
  geom_histogram(binwidth = 1, color = "white") +
  facet_wrap(~site_cluster, scales = "free_y", ncol = 2) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Temporal Distribution for Top 10 Clusters",
    y = "Number of Samples",
    x = "Year"
  ) +
  theme(legend.position = "bottom")

# ==============================================================================
# MODULE 13: CARTOGRAPHIC MASTER MAPS (Dynamic Zoom & Anti-Overplotting)
# ==============================================================================
library(tidyverse)
library(sf)
library(ggspatial) 
library(patchwork) 

cat("\n--- Preparing Map Data ---\n")

# A. Force the dataset into standard Longitude/Latitude (WGS84 / EPSG:4326)
map_data_sf <- oslo_final_sf %>% st_transform(crs = 4326)

# B. Map 1 Data: Extract unique geographic centroids for ALL clusters
df_clusters_unique <- map_data_sf %>%
  group_by(site_cluster) %>%
  summarise(geometry = st_centroid(st_combine(geometry)), .groups = "drop")

# C. Map 2 Data: Extract unique clusters BUT keep the Region label
df_regions_unique <- map_data_sf %>%
  group_by(site_cluster, Fjord_Part) %>%
  summarise(geometry = st_centroid(st_combine(geometry)), .groups = "drop") %>%
  # Strip any accidental hidden spaces in the names that break color mapping
  mutate(Fjord_Part = trimws(as.character(Fjord_Part)))

# D. Dynamic Bounding Box: Calculate the exact limits of your data and add a 5% buffer!
bbox <- st_bbox(df_clusters_unique)
map_boundaries <- coord_sf(
  xlim = c(bbox["xmin"] - 0.05, bbox["xmax"] + 0.05), 
  ylim = c(bbox["ymin"] - 0.05, bbox["ymax"] + 0.05), 
  crs = 4326, 
  expand = FALSE
)

# E. Set custom, high-contrast colors (Matched exactly to trimmed names)
custom_region_colors <- c(
  "Drammensfjord"   = "#e41a1c", 
  "Inner Oslofjord" = "#377eb8", 
  "Outer East"      = "#4daf4a", 
  "Outer West"      = "#984ea3"  
)

# ---------------------------------------------------------
# 1. Map 1: The Site Clusters (Scrambled Confetti Colors)
# ---------------------------------------------------------
# First, figure out exactly how many unique clusters you have
n_clusters <- length(unique(df_clusters_unique$site_cluster))

# Generate a large palette of colors and randomly shuffle their order!
# We use set.seed() so the "random" colors stay exactly the same every time you run the script.
set.seed(42) 
scrambled_palette <- sample(scales::hue_pal()(n_clusters))

p_clusters <- ggplot() +
  
  # Fetch the live basemap
  annotation_map_tile(type = "cartolight", zoom = 10, alpha = 0.8) +
  
  # Plot UNIQUE clusters colored by ID
  geom_sf(data = df_clusters_unique, aes(color = as.factor(site_cluster)), size = 1.8, alpha = 0.9) +
  
  # THE NEW FIX: Apply the randomized color palette
  scale_color_manual(values = scrambled_palette) +
  
  theme_bw() +
  map_boundaries + 
  labs(title = "A) 1km Site Clusters (N=104)") +
  
  # Explicitly turn off the guide/legend for the color aesthetic
  guides(color = "none") + 
  
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_blank(), 
    axis.text = element_text(size = 8)
  )

# ---------------------------------------------------------
# 2. Map 2: The Four Regions Colored
# ---------------------------------------------------------
p_regions <- ggplot() +
  
  # Fetch the live basemap
  annotation_map_tile(type = "cartolight", zoom = 10, alpha = 0.8) +
  
  # Plot the UNIQUE clusters colored by region (Prevents black overplotting)
  geom_sf(data = df_regions_unique, aes(color = Fjord_Part), size = 2.0, alpha = 0.9) +
  
  scale_color_manual(values = custom_region_colors) +
  
  theme_bw() +
  map_boundaries + 
  labs(title = "B) Analytical Fjord Regions") +
  
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_blank(), 
    axis.text = element_text(size = 8),
    legend.position = "right",
    legend.title = element_blank(), 
    legend.text = element_text(size = 9)
  )

# ---------------------------------------------------------
# 3. Adding Cartographic Elements
# ---------------------------------------------------------
cat("\n--- Finalizing Cartography ---\n")

p_clusters_final <- p_clusters +
  annotation_scale(
    location = "bl", 
    width_hint = 0.4, 
    unit_category = "metric", 
    bar_cols = c("black", "white"),
    text_cex = 0.8, 
    pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm")
  ) +
  annotation_north_arrow(
    location = "tr", 
    which_north = "true",
    style = north_arrow_fancy_orienteering(),
    pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
    height = unit(1, "cm"), width = unit(1, "cm")
  )

p_regions_final <- p_regions + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ---------------------------------------------------------
# 4. Combine and Print
# ---------------------------------------------------------
master_map <- p_clusters_final + p_regions_final + plot_layout(ncol = 2)

cat("\nPrinting Final Map...\n")
print(master_map)

