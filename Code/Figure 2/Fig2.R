# Code to create Figure 2 of Manuscript

# Author(s): Charlotte Ward
# Version: YYYY-MM-DD

# Load Pkgs
library(sf)
library(terra)
library(ggplot2)
library(IsoriX)
library(tidyverse)
library(igraph)
library(vegan)
library(ggspatial)
#Un-hash the below to download the glatos R package to remove false detections
# library(remotes)
# install.packages('glatos', repos = c('https://ocean-tracking-network.r-universe.dev', 'https://cloud.r-project.org'))
library(glatos)

file_path <- getwd()

# load data
waterbodies <- st_read(file.path(file_path, "Data/Lakes_Streams/WaterbodyPolygons.shp"))
watercourses <- st_read(file.path(file_path, "Data/Lakes_Streams/WatercourseLines.shp"))

# Ensure CRS is consistent (WGS84)
waterbodies <- st_transform(waterbodies, crs = 4326)
watercourses <- st_transform(watercourses, crs = 4326)

# Load receiver data
metadata <- read_csv(file.path(file_path, "Data/receivers.csv"))

# Convert receiver locations to sf object
receiver_sf <- st_as_sf(
  metadata,
  coords = c("deploy_long", "deploy_lat"),
  crs = 4326
)

# Load detection data
dets <- read_csv(file.path(file_path, "Data/detections_clean_alldata.csv"))

# Remove specified transmitters
transmitters_to_remove <- c("34905", "32999", "42348")
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
GNIPData <- read.csv(file.path(file_path, "Data/baselines_coords_d13C.csv"))
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

# Figure 2: Isoscape + Network + Receiver Detection Frequencies
ggplot() +
  geom_raster(data = masked_isoscape_df, aes(x = x, y = y, fill = mean)) +
  scale_fill_gradientn(
    colors = c("#b34332", "#F06C57", "#eb8d7f", "#b6cedb","#5494b8","#2c6c91"),
    name = expression(δ^13*C)
  ) +
  geom_sf(data = cropped_waterbodies, fill = NA, size = 0.2) +
  geom_sf(
    data = receiver_sf,
    aes(size = size),
    shape = 21,  
    color = "black",
    fill = "white",
    alpha = 0.7,
    stroke = 1.5
  ) +
  scale_size_continuous(
    name = "Receiver Detections",
    range = c(0.5, 4),  
    breaks = scales::rescale(c(0, 100, 500, 5000, 10000, 20000), to = c(0.5, 3)),  
    labels = c("0", "100", "500", "5000", "10000", "20000")  
  ) +
  annotation_scale(location = "bl", width_hint = 0.3, bar_units = "m", unit_category = "metric", height = unit(0.3, "cm")) +
  annotation_north_arrow(location = "br", style = north_arrow_minimal, height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    axis.ticks = element_blank(), 
    axis.text = element_text(size = 12),   
    axis.title = element_text(size = 14, face = "bold"), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10),  
    aspect.ratio = 9/10  
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  )





