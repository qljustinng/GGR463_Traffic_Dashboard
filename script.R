# IMPORTS ----------------------------------------------------------------------
# Spatial Data manipulation
library(sf)
library(raster)
library(sp)
library(rgdal)
library(maptools)

# Data manipulation
library(tidyverse)

# Statistical Tools
library(gstat)
library("SpatialKDE")
library(Metrics)

# Mapping library
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(leaflet)
library(htmlwidgets)

# DATA -------------------------------------------------------------------------
RLC <- st_read("cleaned/red_light_cameras.shp")
KSI <- st_read("cleaned/KSI.shp")
#traffic <- st_read("cleaned/traffic_2000-2029_2.shp")
boundary <- readOGR("cleaned/citygcs.ward_2018_wgs84.shp")

boundary <- spTransform(boundary, "+proj=longlat +datum=WGS84 +no_defs")

# PARAMETERS -------------------------------------------------------------------

# 
# Ascending order
dateSegments <- as.POSIXct(c("2006-01-01", "2007-01-01", "2008-01-01", 
                             "2009-01-01", "2010-01-01", "2011-01-01",
                             "2012-01-01", "2013-01-01", "2014-01-01",
                             "2015-01-01", "2016-01-01", "2017-01-01",
                             "2018-01-01", "2019-01-01", "2020-01-01"))

# IDW
K <- 2 # weight of distance for IDW

# KDE
cell_size <- 500 # larger number == lower resolution
band_width <- 3000 # larger number == larger area of effect

# COMMON CODE ------------------------------------------------------------------



# Separate by time segments ----------------------------------------------------

# Seperate by time segments
dateMap <- function (date) {
  if (is.na(date)) {
    return(-1)
  }
  
  for(j in 1:length(dateSegments)) {
    if (date < dateSegments[j]) {
      return(j)
    }
  }
  return(length(dateSegments) + 1)
}

dateMap2 <- function (dateUnformattedVector) {
  output <- list()
  print(length(dateUnformattedVector))
  for (i in 1:length(dateUnformattedVector)) {
    date <- as.POSIXct(substr(dateUnformattedVector[i], 1, 10)[1])
    output <- append(output, as.integer(dateMap(date)))
  }
  
  return(output)
}

# Collisions
KSI <- KSI %>% mutate(dateClass = as.integer(dateMap2(DATE)))
KSI_unique <- unique(KSI$dateClass)

for (i in 1:length(KSI_unique)) {
  tmp <- KSI %>% filter(dateClass == KSI_unique[i])
  st_write(tmp, paste0("KSI/", "KSI", as.character(KSI_unique[i]), ".shp"))
}

# Traffic Volume
# This section uses QGIS, 602729 observations, more efficient in QGIS

rm(dateMap, dateMap2, KSI, RLC, i)

# Creating IDW interpolated surfaces -------------------------------------------

# Shared between all IDW rasters

# Union of boundary file into single polygon
full_boundary <- unionSpatialPolygons(boundary, IDs=boundary@data$FEATURE_CO)
#plot(full_boundary)

# Create empty grid within extent of boundary
grid <- makegrid(full_boundary, cellsize = .002)
coordinates(grid) <- ~ x1+x2
proj4string(grid) <- full_boundary@proj4string@projargs
#plot(grid)

# only keep points inside grid
grid_in_poly <- sp::over(grid, full_boundary)
boundary_grid <- grid[!is.na(grid_in_poly),]
#plot(boundary_grid)

# get all traffic filenames
filenames <- list.files("TRAFFIC", pattern="*.shp", full.names=TRUE)

RMSE <- list()

# Complete IDW for all traffic files
for(traffic_filename in filenames) {
  traffic <- readOGR(traffic_filename)
  
  traffic_unique <- remove.duplicates(traffic)
  traffic_unique$total_traf <- traffic_unique$total_traf %>% replace_na(0)

  # Create a IDW model using a K value of 2 since it is a default
  traffic_idw <- idw(total_traf ~ 1,
                     traffic_unique,
                     boundary_grid,
                     idp = K
  )
  
  traffic_idw_frame = as.data.frame(traffic_idw)
  
  # Visualize IDW
  # spplot(traffic_idw, "var1.pred")
  
  # Residuals
  LOOCV <- krige.cv(total_traf~1, 
                    traffic_unique, 
                    nfold = nrow(traffic_unique),
                    set = list(idp = K)
  )
  
  # Visualize residuals
  #spplot(LOOCV, "residual")
  
  # Compute RMSE
  RMSE <- c(RMSE, LOOCV@data$residual ^ 2 %>%
    mean() %>%
    sqrt())
  
  # save as shapefile (optional)
  #coordinates(traffic_idw_frame) = ~ x1+x2
  #proj4string(traffic_idw_frame) <- traffic_idw@proj4string@projargs
  #raster::shapefile(traffic_idw_frame, paste0("TRAFFIC_IDW/TRAFFIC_IDW", 
  #                                            substr(traffic_filename, 16, 
  #                                                   nchar(traffic_filename))))
  
  # save as tif
  traffic_idw_rast <- rasterFromXYZ(traffic_idw)
  writeRaster(traffic_idw_rast, paste0("TRAFFIC_IDW/TRAFFIC_IDW", 
                                       substr(traffic_filename, 16, 
                                              nchar(traffic_filename) - 3), 
                                       "tif"))
  
  # Save residuals
  raster::shapefile(LOOCV, paste0("TRAFFIC_IDW_RESIDUALS/TRAFFIC_IDW_RESIDUALS", 
                                  substr(traffic_filename, 16, 
                                         nchar(traffic_filename))))
}

RMSE_values <- data.frame(IDW_filename = filenames,
                 RMSE_value = unlist(RMSE))

write.csv(RMSE_values, "TRAFFIC_IDW_RMSE.csv", row.names=TRUE)

rm(boundary_grid, full_boundary, grid, LOOCV, RMSE, RMSE_values, traffic, 
   traffic_idw, traffic_idw_rast, traffic_unique, filenames, i, K, 
   traffic_filename)

# Creating Weighted Kernel Density ---------------------------------------------

KSI_files <- list.files("KSI", pattern="*.shp", full.names=TRUE)
traffic_idw_files <- list.files("TRAFFIC_IDW", pattern="*.tif", full.names=TRUE)

clean_ksi_filename <- function(name) {
  return(substr(name, nchar("KSI/KSI") + 1, nchar(name) - nchar(".shp")))
}

clean_idw_filename <- function(name) {
  return(substr(name, nchar("TRAFFIC_IDW/TRAFFIC_IDW") + 1, 
                nchar(name) - nchar(".tif")))
}

idw_valid <- map(traffic_idw_files, clean_idw_filename)
KSI_valid <- map(KSI_files, clean_ksi_filename)
valid_time_codes <- intersect(idw_valid, KSI_valid)

# loop over
for(time_code in valid_time_codes) {
  traffic <- raster(paste0("TRAFFIC_IDW/TRAFFIC_IDW", time_code, ".tif"))
  KSI <- readOGR(paste0("KSI/KSI", time_code, ".shp"))
  KSI_transformed <- spTransform(KSI, proj4string(traffic))
  
  # Calculate weight
  KSI_transformed@data$traffic <- 0
  KSI_transformed@data$traffic <- 1 / raster::extract(traffic, KSI)
  KSI_transformed@data$traffic <- KSI_transformed@data$traffic %>% replace_na(0)
  
  # KDE change data to correct format
  KSI_transformed <- spTransform(KSI_transformed, CRS("EPSG:32617"))
  KSI_sf <- as(KSI_transformed, "sf")
  KSI_sf <-  KSI_transformed %>% st_as_sf(coords = c("x", "y"), dim = "XY") %>% 
    select()

  # KDE grid
  kde_grid <- KSI_sf %>%
    #create_grid_rectangular(cell_size = cell_size, side_offset = band_width)
    #create_grid_hexagonal(cell_size = cell_size, side_offset = band_width)
    create_raster(cell_size = cell_size, side_offset = band_width)
  
  # Compute
  kde <- KSI_sf %>%
    kde(band_width = band_width, kernel = "quartic", grid = kde_grid, weights = KSI_transformed@data$traffic)
  
  # Plot rectangular or hexagonal grid
  # tm_shape(c) +
  #   tm_polygons(col = "kde_value", palette = "viridis", title = "KDE") +
  #   tm_shape(KSI_sf) +
  #   tm_bubbles(size = 0.05, col = "red") +
  #   tm_shape(boundary) +
  #   tm_borders("white", lwd = .5)
  
  # Plot raster layer
  # raster_df <- as.data.frame(kde, xy=TRUE)
  # ggplot() +
  #   geom_tile(data = raster_df , aes(x = x, y = y, fill = layer)) + 
  #   scale_fill_distiller(palette = "Spectral") +
  #   geom_sf(data = KSI_sf, colour = "red", fill = NA, size = 0.1)
  
  writeRaster(kde, paste0("KSI_KDE/KSI_KDE", time_code, ".tif"))
}

rm(clean_idw_filename, clean_ksi_filename, band_width, cell_size, KSI_files, 
   traffic_idw_files, time_code, idw_valid, kde, kde_grid, KSI, KSI_sf, 
   KSI_transformed, KSI_valid, raster_df, traffic, valid_time_codes)

# Sample raster values ---------------------------------------------------------

RLC <- st_read("cleaned/red_light_cameras.shp")

# Convert to dataframe (for use in tidyverse)
RLC_df <- as.data.frame(RLC) %>% 
  select(ACTIVATION, geometry)

# clean filename function
clean_kde_filename <- function(filename) {
  return(substr(filename, nchar("KSI_KDE/KSI_KDE") + 1, 
                nchar(filename) - nchar(".tif")))
}

# Get unique time id for KDE and RLC
KDE_files <- list.files("KSI_KDE", pattern="*.tif", full.names=TRUE)
KDE_valid <- KDE_files %>% map(clean_kde_filename)

# Iterate through each raster and extract value to RLC dataframe
for (time_code in KDE_valid) {
  KDE <- raster(paste0("KSI_KDE/KSI_KDE", time_code, ".tif"))
  extract <- raster::extract(KDE, RLC)

  RLC_df <- RLC_df %>%
    add_column(extract) %>%
    rename(!!paste0("tc_", as.character(time_code)) := extract)
}

RLC_sdf <- st_as_sf(RLC_df)

# Plot out data
ggplot() + 
  geom_sf(data = RLC_sdf, aes(color = tc_10)) +
  scale_fill_distiller(palette = "Spectral")

# Save to file
st_write(RLC_sdf, "RLC_HOTSPOT.shp")

rm(clean_kde_filename, extract, KDE_files, time_code, KDE, KDE_valid, RLC, 
   RLC_df, RLC_sdf)

# Linear Regression & Map Creation ---------------------------------------------

# Functions to convert dates to time codes using tidyverse
dateMap <- function (date) {
  if (is.na(date)) {
    return(-1)
  }
  
  for(j in 1:length(dateSegments)) {
    if (date < dateSegments[j]) {
      return(j)
    }
  }
  return(length(dateSegments) + 1)
}

dateMap2 <- function (dateUnformattedVector) {
  output <- list()
  for (i in 1:length(dateUnformattedVector)) {
    date <- as.POSIXct(substr(dateUnformattedVector[i], 1, 10)[1])
    output <- append(output, as.integer(dateMap(date)))
  }
  
  return(output)
}

# Ward boundary
boundary <- readOGR("cleaned/citygcs.ward_2018_wgs84.shp")
boundary <- spTransform(boundary, "+proj=longlat +datum=WGS84 +no_defs")
boundary <- as(boundary, "sf")

# Read in files and add time code
RLC <- st_read("RLC_HOTSPOT.shp")
RLC <- RLC %>% mutate(activeDateClass = as.integer(dateMap2(ACTIVATION)))
RLC <- as.data.frame(RLC)

# Get all column names (time codes)
RLC_cols <- colnames(RLC)
RLC_cols <- RLC_cols[is.na(match(RLC_cols,c("ACTIVATION", "geometry", 
                                            "activeDateClass")))]

# Reformat frames to apply linear regression
RLC_reformat3 <- list()
# find builtin (tidyverse?)
for (i in 1:nrow(RLC)) {
  active_date_class <- as.numeric(RLC[i, "activeDateClass"])
  
  temp_3 <- data.frame(time = NA, accidentRate = NA)[-1, ]
  
  for (col in RLC_cols) {
    time_class_relative <- as.numeric(substr(col, 4, 
                                             nchar(col))) - active_date_class
    temp_3 <- rbind(temp_3, data.frame(time = time_class_relative, 
                                       accidentRate = RLC[i, col]))
  }
  RLC_reformat3[[i]] <- temp_3
}

# Stores linear regression information
multipleRegressionPartialBefore <- list()
multipleRegressionPartialAfter <- list()

# Variables for regression
# y = mx + b, R squared

# R squared
r_squared_before <- list()
r_squared_after <- list()

# m
coeff_before <- list()
coeff_after <- list()

# b
intercept_before <- list()
intercept_after <- list()

# coeff_after - coeff_before
coeff_diff <- list()

# Apply linear regression
for(i in 1:length(RLC_reformat3)) {

  # Two regression lines (before and after)
  # Reformat to dataframe
  df <- map2_dfr(RLC_reformat3[[i]]$time, RLC_reformat3[[i]]$accidentRate, ~ tibble(time = .x, accidentRate = .y))
  
  # Before regression
  # complete regression if possible
  if (length(df[df$time < 0, ]$time) > 1) {
    multipleRegressionPartialBefore[[i]] <- lm(df[df$time < 0, ]$time ~ df[df$time < 0, ]$accidentRate)
    
    if (nrow(summary(multipleRegressionPartialBefore[[i]])$coefficients) == 1) {
      # horizontal line/ not enough data
      coeff_before <- append(coeff_before, -9999)
      intercept_before <- append(intercept_before, -9999)
      r_squared_before <- append(r_squared_before, -9999)
    } else {
      # get slope
      coeff_before <- append(coeff_before, round(summary(multipleRegressionPartialBefore[[i]])$coefficients[2, 1], digits=3))
      intercept_before <- append(intercept_before, round(summary(multipleRegressionPartialBefore[[i]])$coefficients[1, 1], digits=3))
      r_squared_before <- append(r_squared_before, round(summary(multipleRegressionPartialBefore[[i]])$r.squared, digits=3))
    }
  } else {
    # could not do regression
    coeff_before <- append(coeff_before, -9999)
    intercept_before <- append(intercept_before, -9999)
    r_squared_before <- append(r_squared_before, -9999)
  }
  
  # Same as above but for the After regression
  if (length(df[df$time > 0, ]$time) > 1) {
    multipleRegressionPartialAfter[[i]] <- lm(df[df$time > 0, ]$time ~ df[df$time > 0, ]$accidentRate)
    
    if (nrow(summary(multipleRegressionPartialAfter[[i]])$coefficients) == 1) {
      # regression not possible
      coeff_after <- append(coeff_after, -9999)
      intercept_after <- append(intercept_after, -9999)
      r_squared_after <- append(r_squared_after, -9999)
    } else {
      coeff_after <- append(coeff_after, round(summary(multipleRegressionPartialAfter[[i]])$coefficients[2, 1], digits=3))
      intercept_after <- append(intercept_after, round(summary(multipleRegressionPartialAfter[[i]])$coefficients[1, 1], digits=3))
      r_squared_after <- append(r_squared_after, round(summary(multipleRegressionPartialAfter[[i]])$r.squared, digits=3))
    }
  } else {
    coeff_after <- append(coeff_after, -9999)
    intercept_after <- append(intercept_after, -9999)
    r_squared_after <- append(r_squared_after, -9999)
  }
  
  if (coeff_after[[i]] == -9999 | coeff_before[[i]] == -9999) {
    coeff_diff <- append(coeff_diff, -9999)
  } else {
    coeff_diff <- append(coeff_diff, coeff_after[[i]] - coeff_before[[i]]) # negative means improvement
  }
    
}

# Creating layout
layout <- rbind(c(NA,NA,NA,2,2),
                c(NA,NA,NA,2,2),
                c(1,1,1,1,NA),
                c(1,1,1,1,NA),
                c(1,1,1,1,NA),
                c(1,1,1,1,NA))

# Analysis use
for(i in 1:length(RLC_reformat3)) {
  df <- data.frame(long = RLC[["geometry"]][[i]][[1]], lat = RLC[["geometry"]][[i]][[2]]) 
  data_item <- RLC_reformat3[[i]] %>% mutate(beforeAfter = case_when(time < 0 ~ 1, time > 0 ~ 2, time == 0 ~ 0))
  
  grid <- grid.arrange(
    ggplot(data = data_item, aes(x = time, y = accidentRate, group=beforeAfter)) + 
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, aes(color=beforeAfter)) +
      theme(legend.position = "none") +
      labs(title="Accident Rate of Red-light Camera",
           x ="Years since installation", 
           y = "Accident Rate"), 
    
    ggplot() + 
      geom_sf(data = boundary) +
      geom_point(data = df, aes(x = long, y = lat), col = 'blue') +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()),
    
    layout_matrix=layout)
  
  ggsave(paste0("LM/", as.character(i), ".png"), grid)
}

# Mapping use
for(i in 1:length(RLC_reformat3)) {
  data_item <- RLC_reformat3[[i]] %>% mutate(beforeAfter = case_when(time < 0 ~ 1, time > 0 ~ 2, time == 0 ~ 0))
  img <- ggplot(data = data_item, aes(x = time, y = accidentRate, group=beforeAfter)) + 
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, aes(color=beforeAfter)) +
    theme(legend.position = "none") +
    labs(title="Accident Rate of Red-light Camera",
         x ="Years since installation", 
         y = "Accident Rate")
  ggsave(paste0("LM_app/", as.character(i), ".png"), img)
}

# Adding information to map
RLC_spdf <- readOGR("RLC_HOTSPOT.shp")
RLC_spdf@data$coeff_before <- unlist(coeff_before)
RLC_spdf@data$intercept_before <- unlist(intercept_before)
RLC_spdf@data$coeff_after <- unlist(coeff_after)
RLC_spdf@data$intercept_after <- unlist(intercept_after)
RLC_spdf@data$r_squared_before <- unlist(r_squared_before)
RLC_spdf@data$r_squared_after <- unlist(r_squared_after)
RLC_spdf@data$coeff_diff <- unlist(coeff_diff)
RLC_spdf@data$img_no <- unlist(1:nrow(RLC_spdf)) # original IDs

RLC_spdf_filtered <- as.data.frame(RLC_spdf)
RLC_spdf_filtered <- RLC_spdf_filtered %>% mutate(goodData = ifelse(coeff_before == -9999 | 
                                                    intercept_before == -9999 |
                                                    coeff_after == -9999 | 
                                                    intercept_after == -9999, 
                                                  FALSE, TRUE)) %>% filter(goodData == TRUE)

RLC_spdf_filtered <- st_as_sf(RLC_spdf_filtered, coords = c("coords.x1", "coords.x2"), dim = "XY")

# Interactive leaflet map
map <- leaflet() %>%
  addTiles() %>%
  setMaxBounds(lng1 = -80, lat1 = 43.5, lng2=-78.8, lat2=44.1) %>%
  addPolygons(data = boundary, fillOpacity = 0, group = "Ward Boundaries") %>%
  addMarkers(data = RLC_spdf, popup = 
                     paste0("<p> Activiation Date: ", RLC_spdf$ACTIVATION, "<br>",
                            "Equation Before Activation: y = ", RLC_spdf$coeff_before, "x + ", RLC_spdf$intercept_before, "<br>",
                            "R squared Before Activation: ", RLC_spdf$r_squared_before, "<br>",
                            "Equation After Activation: y = ", RLC_spdf$coeff_after, "x + ", RLC_spdf$intercept_after, "<br>",
                            "R squared After Activation: ", RLC_spdf$r_squared_after, "<br>",
                            
                            "<a href='LM_app/", as.character(c(1:length(RLC_reformat3))), ".png'>",
                            "<img width='250px' src='LM_app/", as.character(c(1:length(RLC_reformat3))), ".png' />",
                            "</a>"), group = "All Red-light Cameras") %>%
  
  addMarkers(data = RLC_spdf_filtered, popup = 
               paste0("<p> Activiation Date: ", RLC_spdf_filtered$ACTIVATION, "<br>",
                      "Equation Before Activation: y = ", RLC_spdf_filtered$coeff_before, "x + ", RLC_spdf_filtered$intercept_before, "<br>",
                      "R squared Before Activation: ", RLC_spdf_filtered$r_squared_before, "<br>",
                      "Equation After Activation: y = ", RLC_spdf_filtered$coeff_after, "x + ", RLC_spdf_filtered$intercept_after, "<br>",
                      "R squared After Activation: ", RLC_spdf_filtered$r_squared_after, "<br>",
                      
                      "<a href='LM_app/", RLC_spdf_filtered$img_no, ".png'>",
                      "<img width='250px' src='LM_app/", RLC_spdf_filtered$img_no, ".png' />",
                      "</a>"), group = "Possible Regression") %>%
  
  addLayersControl(
    overlayGroups = c("All Red-light Cameras", "Possible Regression", "Ward Boundaries"),
    options = layersControlOptions(collapsed = FALSE)) %>%
  htmlwidgets::onRender("
        function() {
            $('.leaflet-control-layers-overlays').prepend('<label style=\"text-align:center\"><b>Toronto Red-light cameras</b></label>')
            .append('<label> <hr>Click on each point to see more</label> <label>information. Click on the graph in </label> <label>the popup to enlarge the image.<hr>Quoc Luan Nguyen | 1005855184</label>');
        }
    ") %>%
  setView(-79.358145, 43.731887, zoom=11)

# View (not fully working) and export
map
saveWidget(map, file="app.html")

# Add custom tab icon and name
html_file <- readLines("app.html")
html_file  <- gsub(pattern = "<title>leaflet</title>", replace = "<title>Red-light Camera Traffic Accidents</title><link rel=\"Page icon\" href=\"red_light.ico\">", x = html_file)
writeLines(html_file, con="RLC_app.html")

# Regression information

# create two variable dataframe for histogram
RLC_summary_list <- c(RLC_spdf_filtered$coeff_before, RLC_spdf_filtered$coeff_after, RLC_spdf_filtered$coeff_diff)
RLC_summary_desig <- c(unlist(list(rep("Before", length(RLC_spdf_filtered$coeff_before)))),
                       unlist(list(rep("After", length(RLC_spdf_filtered$coeff_after)))),
                       unlist(list(rep("Difference", length(RLC_spdf_filtered$coeff_before)))))

RLC_summary <- data.frame(data_type = RLC_summary_desig, 
                          value = RLC_summary_list)

# histogram + mean for each coefficient
data_text <- data.frame(
  label = c(paste0("Mean: ", as.character(round(mean(RLC_spdf_filtered$coeff_before), 3))), 
            paste0("Mean: ", as.character(round(mean(RLC_spdf_filtered$coeff_after), 3))), 
            paste0("Mean: ", as.character(round(mean(RLC_spdf_filtered$coeff_diff), 3))),
            paste0("Mean R^2: ", as.character(round(mean(RLC_spdf_filtered$r_squared_before), 3))),
            paste0("Mean R^2: ", as.character(round(mean(RLC_spdf_filtered$r_squared_after), 3)))
            ),
  data_type   = c("Before", "After", "Difference", "Before", "After"),
  value = c(-2000, -2000, -2000, 2000, 2000),
  count = c(45, 45, 45, 45, 45)
)

# output result
img <- ggplot(data = RLC_summary, aes(x = value)) +
              geom_histogram(binwidth = 50, aes(fill=data_type)) + 
              facet_wrap(~factor(data_type, levels=c("Before", "After", "Difference")), ncol = 1) +
              geom_text(data = data_text, mapping = aes(x = value, y = count, label = label)) +
              theme(legend.position = "none")

ggsave("regression.png", img)

RLC_summary_list <- c(RLC_spdf_filtered$coeff_before, RLC_spdf_filtered$coeff_after, RLC_spdf_filtered$coeff_diff)
RLC_summary_desig <- c(unlist(list(rep("Before", length(RLC_spdf_filtered$coeff_before)))),
                       unlist(list(rep("After", length(RLC_spdf_filtered$coeff_after)))),
                       unlist(list(rep("Difference", length(RLC_spdf_filtered$coeff_before)))))

rm(boundary, img, map, multipleRegression, r_squared, RLC, RLC_reformat3, 
   temp_3, x_const, x_const_p, x_int, x_int_p, html_file, dateMap, dateMap2, 
   active_date_class, col, dateSegments, i, RLC_cols, time_class_relative, df, 
   coeff_before, coeff_after, coeff_diff, data_item, data_text, grid, 
   intercept_before, intercept_after, layout, multipleRegressionPartialBefore, 
   multipleRegressionPartialAfter, RLC_spdf, RLC_spdf_filtered, 
   r_squared_before, r_squared_after, RLC_summary, RLC_summary_desig, 
   RLC_summary_list)
