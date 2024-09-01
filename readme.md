 
# APRIL Spatial Database Project

This project focuses on the analysis of spatial data, specifically land use and elevation in Indonesia. The project includes Python scripts for data loading and report generation, an R script for spatial analysis, and a SQL file for querying linked spatial data.

## Prerequisites

Before running this project, ensure that you have the following installed:
- Docker: For setting up services via `docker-compose`.
- Python 3.8+: For running data processing and report generation scripts.
- R 4.0+: For executing the spatial analysis R script.
- Access to Google Cloud Platform: For authentication using the provided service account credentials.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/APRIL_spatial_db.git
   cd APRIL_spatial_db
   ```

2. Set up Docker services:
   ```bash
   docker-compose up
   ```

3. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Install R dependencies:
   You can install required R packages by running:
   ```R
   install.packages(c("sf", "sp", "raster", "ggplot2"))
   ```

## Code Structure and Explanation

### 1. `scripts/python/load_data.py`

This Python script is responsible for loading and processing the raw spatial data stored in CSV and GeoPackage formats. 

- **Purpose**: To load data from files like `lu.csv` and `lu.gpkg` into a structured format suitable for analysis.
- **Key Components**:
  - Reads data from the `data/` directory.
  - Cleans and formats the data for use in downstream analyses.
  - Likely interacts with a spatial database or processes GeoPackage files using Python libraries like `geopandas`.
  
  Example of how to run:
  ```bash
  python scripts/python/load_data.py
  ```

### 2. `scripts/python/load_report.py`

This Python script generates reports based on the processed data. It outputs the results in CSV format and can create visualizations.

- **Purpose**: To produce a summary report and various visual outputs such as maps and histograms.
- **Key Components**:
  - Utilizes libraries such as `pandas`, `matplotlib`, and `seaborn` to generate visualizations.
  - Generates summary reports, including statistical summaries and plots based on the spatial data.
  
  Example of how to run:
  ```bash
  python scripts/python/load_report.py
  ```

### 3. `scripts/R/spatial_analysis_indonesia.r`

This R script performs spatial analysis on Indonesian land use data, elevation, and other geographic metrics.

- **Purpose**: To carry out spatial data analysis and generate high-resolution maps and visualizations.
- **Key Components**:
  - Reads GeoPackage and other spatial data formats.
  - Uses R libraries like `sf`, `ggplot2`, and `raster` to manipulate and visualize geographic data.
  - Outputs files such as `geometries_map.png` and `high_res_mean_elevation_map.png` stored in the `results/` directory.
  
  Example of how to run:
  ```bash
  Rscript scripts/R/spatial_analysis_indonesia.r
  ```

### 4. `query/view_linked_data.sql`

This SQL script defines a query to view linked spatial data in the database.

- **Purpose**: To retrieve and visualize linked spatial data.
- **Key Components**:
  - Likely used with a PostgreSQL/PostGIS database.
  - Fetches data related to spatial geometries, and more for further analysis.
  
  Example of how to run:
  ```sql
  psql -d your_database -f query/view_linked_data.sql
  ```

## Data

The project uses the following data files located in the `data/` directory:
- **`lu.csv`**: A CSV file containing land use data.
- **`lu.gpkg`**: A GeoPackage file with spatial data.

These datasets are essential for performing spatial analysis and generating reports.

## Results

The results of the analysis are stored in the `results/` directory, including:
- **`geometries_map.png`**: A map showing spatial geometries.
- **`high_res_mean_elevation_map.png`**: A high-resolution elevation map.
- **`mean_elevation_histogram.png`**: A histogram of mean elevations.
- **`[YEAR]_ndvi_time_series.png`**: A time series plot of NDVI.
- **`summary_report.csv`**: A CSV file summarizing the findings.

## Usage

To execute the full workflow, follow these steps:

1. **Load and process the data**:
   ```bash
   python scripts/python/load_data.py
   ```

2. **Generate reports and visualizations**:
   ```bash
   python scripts/python/load_report.py
   ```

3. **Perform spatial analysis in R**:
   ```bash
   Rscript scripts/R/spatial_analysis_indonesia.r
   ```

4. **Query linked spatial data**:
   Use the SQL script to view data from your spatial database:
   ```bash
   psql -d your_database -f query/view_linked_data.sql
   ```
 