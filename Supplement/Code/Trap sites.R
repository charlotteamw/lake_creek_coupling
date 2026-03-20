# Code to create Trap Site Map (Figure S1)
# Author(s): Charlotte Ward

library(sf)
library(ggplot2)
library(tidyverse)
library(ggspatial)

file_path <- getwd()

# ---------------------------------------------------------------------------
# Load basemap layers
# ---------------------------------------------------------------------------
waterbodies  <- st_read(file.path(file_path, "Data/Lakes_Streams/WaterbodyPolygons.shp"))
watercourses <- st_read(file.path(file_path, "Data/Lakes_Streams/WatercourseLines.shp"))

waterbodies  <- st_transform(waterbodies,  crs = 4326)
watercourses <- st_transform(watercourses, crs = 4326)

# ---------------------------------------------------------------------------
# Parse the wpt block out of the Garmin multi-section CSV
# ---------------------------------------------------------------------------
raw_lines <- readLines(file.path(file_path, "Data/trap_sites.csv"))

wpt_start  <- which(raw_lines == "wpt") + 1   # header row is one after "wpt"
# next blank line after wpt block marks the end
wpt_end    <- which(raw_lines[(wpt_start + 1):length(raw_lines)] == "")[1] + wpt_start - 1

wpt_block  <- raw_lines[wpt_start:wpt_end]
traps_raw  <- read_csv(paste(wpt_block, collapse = "\n"), show_col_types = FALSE)

traps <- traps_raw %>%
  select(lat, lon, name) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(name)) %>%
  mutate(
    trap_type = case_when(
      str_starts(name, "CR")          ~ "Minnow Trap (Creek)",
      str_starts(name, "LK")          ~ "Minnow Trap (Lake)",
      str_starts(name, "Directional") ~ "Trap Net",
      TRUE                            ~ "Other"
    )
  )

traps_sf <- st_as_sf(traps, coords = c("lon", "lat"), crs = 4326)

# ---------------------------------------------------------------------------
# Crop basemap to trap extent + small buffer
# ---------------------------------------------------------------------------
map_bbox <- st_bbox(traps_sf)
buffer   <- 0.005
xmin <- as.numeric(map_bbox["xmin"]) - buffer
xmax <- as.numeric(map_bbox["xmax"]) + buffer
ymin <- as.numeric(map_bbox["ymin"]) - buffer
ymax <- as.numeric(map_bbox["ymax"]) + buffer

crop_extent <- st_polygon(list(rbind(
  c(xmin, ymin),
  c(xmax, ymin),
  c(xmax, ymax),
  c(xmin, ymax),
  c(xmin, ymin)
))) %>%
  st_sfc(crs = 4326)

cropped_waterbodies  <- st_intersection(waterbodies,  crop_extent)
cropped_watercourses <- st_intersection(watercourses, crop_extent)

# ---------------------------------------------------------------------------
# Colour + shape palette
# ---------------------------------------------------------------------------
trap_colours <- c(
  "Minnow Trap (Creek)" = "#F06C57",
  "Minnow Trap (Lake)"  = "#4A90B8",
  "Trap Net"            = "#636262"
)

trap_shapes <- c(
  "Minnow Trap (Creek)" = 21,
  "Minnow Trap (Lake)"  = 21,
  "Trap Net"            = 24
)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
ggplot() +
  geom_sf(data = cropped_waterbodies,  fill = "#c6dcec", colour = "grey40", linewidth = 0.3) +
  geom_sf(data = cropped_watercourses, colour = "#7aafd4", linewidth = 0.5) +
  geom_sf(
    data   = traps_sf,
    aes(fill = trap_type, shape = trap_type),
    colour = "black",
    size   = 3,
    stroke = 0.8,
    alpha  = 0.9
  ) +
  geom_label_repel(
    data         = traps_sf,
    aes(label = name, geometry = geometry),
    stat         = "sf_coordinates",
    size         = 2.5,
    nudge_y      = 0.001,
    box.padding  = 0.3,
    point.padding = 0.2,
    segment.colour = "grey50",
    segment.size   = 0.3,
    min.segment.length = 0.1,
    colour       = "black",
    fill         = alpha("white", 0.7),
    label.size   = NA
  ) +
  scale_fill_manual(values  = trap_colours, name = "Trap Type") +
  scale_shape_manual(values = trap_shapes,  name = "Trap Type") +
  annotation_scale(
    location      = "bl",
    width_hint    = 0.3,
    unit_category = "metric",
    height        = unit(0.3, "cm")
  ) +
  annotation_north_arrow(
    location = "br",
    style    = north_arrow_minimal,
    height   = unit(0.5, "cm"),
    width    = unit(0.5, "cm")
  ) +
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.ticks   = element_blank(),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    aspect.ratio = 9/10
  ) +
  labs(x = "Longitude", y = "Latitude")

ggsave("FigureS1_TrapSites.png", width = 10, height = 9, dpi = 600)
