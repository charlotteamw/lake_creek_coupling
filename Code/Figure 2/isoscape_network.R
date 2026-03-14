# Load required libraries
library(sf)
library(terra)
library(ggplot2)
library(IsoriX)
library(glatos)
library(tidyverse)
library(igraph)
library(vegan)
library(ggspatial)

# Load and transform shapefiles
waterbodies <- st_read("/Users/charlotteward/Documents/algonquin_minnow/Final_Figures/final_figures_data/Lakes_Streams 2/WaterbodyPolygons.shp")
watercourses <- st_read("/Users/charlotteward/Documents/algonquin_minnow/Final_Figures/final_figures_data/Lakes_Streams 2/WatercourseLines.shp")

# Ensure CRS is consistent (WGS84)
waterbodies <- st_transform(waterbodies, crs = 4326)
watercourses <- st_transform(watercourses, crs = 4326)

# Load receiver data
metadata <- read_csv("/Users/charlotteward/Documents/algonquin_minnow/Final_Figures/final_figures_data/receivers.csv")

# Convert receiver locations to sf object
receiver_sf <- st_as_sf(
  metadata,
  coords = c("deploy_long", "deploy_lat"),
  crs = 4326
)

# Load detection data
dets <- read_csv("/Users/charlotteward/Documents/algonquin_minnow/Telemetry Data Processing/detections_clean_alldata.csv")

# Remove specified transmitters
transmitters_to_remove <- c("1594949", "1576137")
detections_filtered <- dets %>%
  filter(!transmitter_id %in% transmitters_to_remove) %>%
  filter(!is.na(receiver_sn))

# Compute receiver detection frequency
receiver_freq <- detections_filtered %>%
  group_by(receiver_sn) %>%
  summarise(freq = n())

# Merge detection frequency with receiver locations
receiver_sf <- receiver_sf %>%
  left_join(receiver_freq, by = c("receiver_sn" = "receiver_sn"))

# Replace NA frequencies with 0 for consistency
receiver_sf$freq[is.na(receiver_sf$freq)] <- 0

# Scale receiver sizes by detection frequency
receiver_sf$size <- scales::rescale(receiver_sf$freq, to = c(2, 7)) # Adjust size range

# Generate network edges
make_network_data <- function(dat) {
  ids <- unique(dat$transmitter_id)
  
  from <- character()
  to <- character()
  fish <- character()
  
  for (i in seq_along(ids)) {
    subs <- dat[dat$transmitter_id == ids[i], ]
    if (nrow(subs) < 2) next
    
    rle_values <- rle(subs$receiver_sn)$values
    from <- append(from, rle_values[-length(rle_values)])
    to <- append(to, rle_values[-1])
    fish <- append(fish, rep(ids[i], length(rle_values) - 1))
  }
  
  individual.moves <- data.frame(from, to, fish, stringsAsFactors = FALSE)
  moves.matrix <- table(individual.moves[, 1:2])
  moves <- as.data.frame(moves.matrix)
  moves <- moves[moves$Freq != 0, ]
  
  return(list(moves = moves, moves.matrix = moves.matrix, individual.moves = individual.moves))
}

full.network <- make_network_data(dat = detections_filtered)

# Add edges as spatial lines
edges <- full.network$moves %>%
  mutate(from = as.character(from), to = as.character(to)) %>%
  rowwise() %>%
  mutate(geometry = st_sfc(st_linestring(
    rbind(
      st_coordinates(receiver_sf[receiver_sf$receiver_sn == from, "geometry"])[1, ],
      st_coordinates(receiver_sf[receiver_sf$receiver_sn == to, "geometry"])[1, ]
    )))) %>%
  st_sf(crs = 4326)

# Load and process the isoscape
GNIPData <- read.csv("/Users/charlotteward/Documents/algonquin_minnow/Final_Figures/final_figures_data/baselines_coords_d13C.csv")
GNIPData$source_ID <- factor(paste("site", GNIPData$lat, GNIPData$long, GNIPData$elev, sep = "_"))

GNIPData_agg <- prepsources(data = GNIPData,
                            long_min = -78.50, long_max = -78.30,
                            lat_min = 45.50, lat_max = 45.85)

# Fit isotopic model
ModelFit <- isofit(data = GNIPData_agg,
                   mean_model_fix = list(elev = FALSE, lat_abs = TRUE))

# Define raster extent
new_extent <- c(xmin = -78.50, xmax = -78.30, ymin = 45.50, ymax = 45.85)

# Create raster template with correct resolution
raster_template <- rast(
  xmin = new_extent["xmin"],
  xmax = new_extent["xmax"],
  ymin = new_extent["ymin"],
  ymax = new_extent["ymax"],
  resolution = 0.00005,  
  crs = "+proj=longlat +datum=WGS84"
)

# Rasterize and mask to waterbodies
waterbodies_raster <- rasterize(vect(waterbodies), raster_template, field = 1, background = NA)
clipped_raster <- mask(waterbodies_raster, vect(waterbodies))

# Create structural raster for isoscape modeling
structural_raster <- prepraster(
  raster = clipped_raster,
  isofit = ModelFit,
  aggregation_factor = 1
)

# Build the isoscape
MinnowIsoscape <- isoscape(
  raster = structural_raster,
  isofit = ModelFit
)

# Mask the isoscape to match **only** the waterbody areas
masked_isoscape <- mask(MinnowIsoscape$isoscapes$mean, vect(waterbodies))

# Convert masked isoscape to a dataframe for ggplot2
masked_isoscape_df <- as.data.frame(masked_isoscape, xy = TRUE, na.rm = TRUE)

# Crop waterbodies and watercourses to match the raster extent
raster_extent <- st_as_sf(as.polygons(ext(masked_isoscape), crs = st_crs(waterbodies)$wkt))
cropped_waterbodies <- st_intersection(waterbodies, raster_extent)
cropped_watercourses <- st_intersection(watercourses, raster_extent)
# Define a new column for receiver size mapping (ensures correct legend)
receiver_sf <- receiver_sf %>%
  mutate(size = scales::rescale(freq, to = c(0.5, 3)))  # Adjust min/max scaling

# Final Plot: Isoscape + Network + Receiver Detection Frequencies
ggplot() +
  # Add the isoscape raster layer (properly masked)
  geom_raster(data = masked_isoscape_df, aes(x = x, y = y, fill = mean)) +
  scale_fill_gradientn(
    colors = c("#2166ac", "#67a9cf", "#d1e5f0", "#dea483", "#cf6b4c", "#b2182b"),
    name = expression(δ^13*C)
  ) +
  # Add waterbody boundaries
  geom_sf(data = cropped_waterbodies, fill = NA, size = 0.2) +
  # Add receiver locations with correct size scaling
  geom_sf(
    data = receiver_sf,
    aes(size = size),
    shape = 21,  # Hollow circle
    color = "black",
    fill = "white",
    alpha = 0.7,
    stroke = 1.5
  ) +
  # Correct legend for receiver detections
  scale_size_continuous(
    name = "Receiver Detections",
    range = c(0.5, 4),  # Ensures correct min/max size scaling
    breaks = scales::rescale(c(0, 100, 500, 5000, 10000, 20000), to = c(0.5, 3)),  # Map real detection values
    labels = c("0", "100", "500", "5000", "10000", "20000")  # Legend matches real detection counts
  ) +
  annotation_scale(location = "bl", width_hint = 0.3, bar_units = "m", unit_category = "metric", height = unit(0.3, "cm")) +
  annotation_north_arrow(location = "br", style = north_arrow_minimal, height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),   # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text = element_text(size = 12),   # Increase axis text size
    axis.title = element_text(size = 14, face = "bold"), # Increase axis label size
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Centered, larger title
    legend.title = element_text(size = 12), # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    aspect.ratio = 9/10  # Keeps proportions correct
  ) +
  # Add axis labels
  labs(
    x = "Longitude",
    y = "Latitude"
  )
