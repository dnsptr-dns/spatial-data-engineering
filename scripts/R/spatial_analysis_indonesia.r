# Install necessary packages if not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# Install and load required libraries
install_if_missing("geodata")
install_if_missing("terra")
install_if_missing("sf")
install_if_missing("DBI")
install_if_missing("RPostgres")
install_if_missing("ggplot2")
install_if_missing("grid")
install_if_missing("ggrepel")

# Function to download, process, and upload spatial data for a given country
process_spatial_data <- function(country_name = "Indonesia", 
                                 db_host = "localhost", 
                                 db_port = 5434, 
                                 db_name = "april_spatial_db", 
                                 db_user = "april", 
                                 db_password = "2030april") {
  # Step 1: Download the Spatial Vector Database for the specified country
  vector_data <- gadm(country = country_name, level = 1, path = tempdir())
  
  # Step 2: Download the Spatial Raster Database (elevation) for the specified country
  raster_data <- elevation_30s(country = country_name, path = tempdir())

  # Step 3: Check and Align CRS if Needed
  print(crs(vector_data))
  print(crs(raster_data))

  # If CRS are different, reproject vector data to match raster data CRS
  if (!identical(crs(vector_data), crs(raster_data))) {
    vector_data <- terra::project(vector_data, crs(raster_data))
    cat("CRS of vector data reprojected to match the raster data.\n")
  }

  # Convert vector data to sf format for ggplot2
  vector_data_sf <- st_as_sf(vector_data)

  # Reproject to a suitable projected coordinate system (e.g., UTM Zone for Indonesia)
  projected_crs <- 32750  # UTM Zone 50S for Indonesia, adjust if necessary
  vector_data_sf <- st_transform(vector_data_sf, crs = projected_crs)

  # Step 4: Extract Statistics from the Raster to the Vector
  extracted_stats <- terra::extract(raster_data, vector_data, fun = mean, na.rm = TRUE)
  print(extracted_stats)

  # Correctly extract and clean the mean elevation data
  mean_elevation <- as.numeric(extracted_stats$IDN_elv_msk)
  mean_elevation <- mean_elevation[!is.na(mean_elevation)]
  vector_data_sf$mean_elevation <- mean_elevation

  # Step 5: Create and Save a Histogram of Mean Elevation Values
  if (length(mean_elevation) > 0) {
    hist_data <- data.frame(mean_elevation = mean_elevation, region = vector_data_sf$NAME_1)
    
    # Calculate the frequency distribution
    hist_data$bin <- cut(hist_data$mean_elevation, breaks = seq(min(mean_elevation), max(mean_elevation), by = 20), right = FALSE)
    freq_distribution <- as.data.frame(table(hist_data$bin))
    colnames(freq_distribution) <- c("mean_elevation_bin", "frequency")
     
    
    # Create the histogram plot
    p <- ggplot(hist_data, aes(x = mean_elevation)) +
      geom_histogram(binwidth = 20, fill = "skyblue", color = "green") +
      labs(title = paste("Histogram of Mean Elevation in", country_name),
          x = "Mean Elevation (meters)", y = "Frequency") +
      theme_minimal() +
      theme(plot.title = element_text(size = 16, face = "bold", color = "white"),
            axis.title = element_text(size = 14, color = "white"),
            axis.text = element_text(size = 12, color = "white"),
            plot.margin = margin(10, 10, 30, 10))
    
    # Save the histogram as a PNG file
    ggsave("results/mean_elevation_histogram.png", plot = p, width = 10, height = 8, units = "in")
    
    # Display inference text on the plot
    inference_text <- "The histogram shows the distribution of mean elevation values across Indonesia's administrative regions."
    grid.text(inference_text, x = 0.5, y = 0.03, gp = gpar(fontsize = 12, col = "white"), draw = TRUE)
    
  } else {
    cat("No valid elevation data available for plotting.\n")
  }



  # Step 6: Create a High-Resolution Map of Mean Elevation Values
  map_plot <- ggplot(data = vector_data_sf) +
    geom_sf(aes(fill = mean_elevation), color = "gray70") +
    geom_text_repel(aes(label = NAME_1, geometry = geometry),
                    stat = "sf_coordinates", size = 2.5, color = "white",
                    max.overlaps = Inf, force = 5, point.padding = 0.5, min.segment.length = 0) +
    scale_fill_viridis_c(option = "viridis", name = "Mean Elevation (meters)",
                         limits = c(min(vector_data_sf$mean_elevation, na.rm = TRUE),
                                    max(vector_data_sf$mean_elevation, na.rm = TRUE))) +
    labs(title = paste("Mean Elevation across Administrative Regions in", country_name),
         subtitle = "Distribution of mean elevation values",
         caption = "Source: Spatial data from geodata and elevation data") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_text(size = 12, face = "bold", color = "white"),
          legend.text = element_text(size = 10, color = "white"),
          plot.title = element_text(size = 16, face = "bold", color = "white"),
          plot.subtitle = element_text(size = 14, color = "white"),
          plot.caption = element_text(size = 10, color = "white"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank())

  # Save the map with higher resolution and larger size
  ggsave("results/high_res_mean_elevation_map.png", plot = map_plot, dpi = 300, width = 14, height = 10, units = "in")
  cat("High-resolution map of mean elevation values created and saved successfully.\n")

  # Step 7: Upload the Vector Database to a Spatial Relational Database
  con <- dbConnect(RPostgres::Postgres(), 
                   host = db_host, port = db_port, dbname = db_name, 
                   user = db_user, password = db_password)

  # Specify the schema and table name
  schema_name <- "staging"
  table_name <- "indonesia_elevation"

  # Upload the vector data to the specified schema and table in the database
  st_write(vector_data_sf, con, Id(schema = schema_name, table = table_name), append = FALSE)

  # Close the database connection
  dbDisconnect(con)
  cat("Spatial vector and raster data have been processed and uploaded to the spatial relational database successfully.\n")
}

# Execute the function with default parameters
process_spatial_data()
