import psycopg2
import geopandas as gpd
import pandas as pd
import ee
import matplotlib.pyplot as plt
from shapely import wkt
from ee import ServiceAccountCredentials
from sqlalchemy import create_engine
from shapely import geometry
import json
import os
import numpy as np
import contextily as ctx

# Initialize Earth Engine with the service account credentials
def initialize_earth_engine(credentials_path, service_account_email):
    """
    Initialize Google Earth Engine with the provided credentials.
    :param credentials_path: Path to the service account credentials JSON file.
    :param service_account_email: Service account email for Earth Engine authentication.
    """
    credentials = ee.ServiceAccountCredentials(service_account_email, credentials_path)
    ee.Initialize(credentials)

# Database connection parameters
def create_db_engine(db_params):
    """
    Create a SQLAlchemy engine for PostgreSQL using the provided parameters.
    :param db_params: Dictionary containing database connection parameters.
    :return: SQLAlchemy engine object.
    """
    return create_engine(f'postgresql://{db_params["user"]}:{db_params["password"]}@{db_params["host"]}:{db_params["port"]}/{db_params["dbname"]}')

# Define the SQL query to select id, keterangan, and their geometries from the database
def fetch_geometries(query, engine):
    """
    Fetch geometries from the database using the provided SQL query.
    :param query: SQL query string to execute.
    :param engine: SQLAlchemy engine object.
    :return: GeoDataFrame containing the fetched geometries.
    """
    return gpd.read_postgis(query, engine, geom_col='geom')

# Function to convert PostGIS geometries to Google Earth Engine compatible geometries
def convert_geom_to_gee(geom):
    """
    Convert a PostGIS Polygon or MultiPolygon geometry to Google Earth Engine Geometry.
    :param geom: Shapely geometry object.
    :return: Earth Engine Geometry object.
    """
    if geom.geom_type == 'Polygon':
        return ee.Geometry.Polygon([list(geom.exterior.coords)])
    elif geom.geom_type == 'MultiPolygon':
        polygons = [list(p.exterior.coords) for p in geom.geoms]
        return ee.Geometry.MultiPolygon(polygons)
    else:
        raise ValueError(f"Unsupported geometry type: {geom.geom_type}. Only Polygon and MultiPolygon are supported.")

# Function to calculate NDVI using Landsat 8 data for a given geometry
def calculate_monthly_ndvi(geometry, year):
    """
    Calculate monthly NDVI for a given geometry using Landsat 8 images from the specified year.
    :param geometry: GEE geometry object.
    :param year: Year for NDVI calculation.
    :return: A dictionary with monthly NDVI values.
    """
    monthly_ndvi = {}
    for month in range(1, 13):
        start_date = f'{year}-{month:02d}-01'
        end_date = f'{year}-{month:02d}-28' if month == 2 else f'{year}-{month:02d}-30'
        
        landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
            .filterDate(start_date, end_date) \
            .filterBounds(geometry) \
            .map(lambda img: img.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI'))
        
        ndvi_image = landsat.median().select('NDVI')
        
        ndvi_value = ndvi_image.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geometry,
            scale=30,
            maxPixels=1e15,
            bestEffort=True
        ).getInfo()
        
        monthly_ndvi[month] = ndvi_value.get('NDVI', None)
    
    return monthly_ndvi

# Function to plot the NDVI time series
def plot_ndvi_time_series(ndvi_results, output_path, year):
    """
    Plot NDVI time series from monthly NDVI values and save to file.
    :param ndvi_results: List of dictionaries containing NDVI values for each category.
    :param output_path: Path to save the plot image.
    """
    plt.figure(figsize=(10, 6))
    for result in ndvi_results:
        label = f"{result['keterangan']}"
        monthly_ndvi = result['monthly_ndvi']
        
        months = list(monthly_ndvi.keys())
        ndvi_values = [monthly_ndvi[month] for month in months]
        
        plt.plot(months, ndvi_values, marker='o', label=label)
    
    plt.xlabel('Month')
    plt.ylabel('NDVI')
    plt.title(f'Monthly NDVI Time Series for Different Categories ({year})')
    plt.xticks(months, [f'{month:02d}' for month in months])  # Format month numbers
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def calculate_annual_ndvi(geometry, start_year, end_year):
    """
    Calculate monthly NDVI for a given geometry using Landsat 8 images from the specified years.
    :param geometry: GEE geometry object.
    :param start_year: Start year for NDVI calculation.
    :param end_year: End year for NDVI calculation.
    :return: A dictionary with monthly NDVI values.
    """
    annual_ndvi = {}
    
    for year in range(int(start_year), int(end_year) + 1):
        for month in range(1, 13):
            start_date = f'{year}-{month:02d}-01'
            end_date = f'{year}-{month:02d}-28' if month == 2 else f'{year}-{month:02d}-30'
            
            # Filter the ImageCollection for the given month and geometry
            landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
                .filterDate(start_date, end_date) \
                .filterBounds(geometry)
            
            # Check if the ImageCollection is empty
            image_count = landsat.size().getInfo()
            if image_count == 0:
                print(f"No images found for {year}-{month:02d}")
                annual_ndvi[f'{year}-{month:02d}'] = None
                continue
            
            # Get a sample image to inspect band names
            sample_image = landsat.first()
            
            # Check if the sample image has the expected bands
            bands = sample_image.bandNames().getInfo()
            if 'SR_B5' not in bands or 'SR_B4' not in bands:
                print(f"Missing bands in collection for {year}-{month:02d}: SR_B5 or SR_B4")
                annual_ndvi[f'{year}-{month:02d}'] = None
                continue
            
            # Calculate NDVI
            landsat_ndvi = landsat.map(lambda img: img.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI'))
            ndvi_image = landsat_ndvi.median().select('NDVI')
            
            ndvi_value = ndvi_image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=geometry,
                scale=30,
                maxPixels=1e15,
                bestEffort=True
            ).getInfo()
            
            annual_ndvi[f'{year}-{month:02d}'] = ndvi_value.get('NDVI', None)
    
    return annual_ndvi


# def calculate_annual_ndvi(geometry, start_year, end_year):
#     """
#     Calculate annual NDVI for a given geometry using Landsat 8 images from the specified years.
#     :param geometry: GEE geometry object.
#     :param start_year: Start year for NDVI calculation.
#     :param end_year: End year for NDVI calculation.
#     :return: A dictionary with annual NDVI values.
#     """
#     annual_ndvi = {}
    
#     for year in range(int(start_year), int(end_year) + 1):
#         start_date = f'{year}-01-01'
#         end_date = f'{year}-12-31'
        
#         # Filter the ImageCollection for the given year and geometry
#         landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
#             .filterDate(start_date, end_date) \
#             .filterBounds(geometry)
        
#         # Check if the ImageCollection is empty
#         image_count = landsat.size().getInfo()
#         if image_count == 0:
#             print(f"No images found for {year}")
#             annual_ndvi[year] = None
#             continue
        
#         # Get a sample image to inspect band names
#         sample_image = landsat.first()
        
#         # Check if the sample image has the expected bands
#         bands = sample_image.bandNames().getInfo()
#         if 'SR_B5' not in bands or 'SR_B4' not in bands:
#             print(f"Missing bands in collection for {year}: SR_B5 or SR_B4")
#             annual_ndvi[year] = None
#             continue
        
#         # Calculate NDVI
#         landsat_ndvi = landsat.map(lambda img: img.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI'))
#         ndvi_image = landsat_ndvi.median().select('NDVI')
        
#         ndvi_value = ndvi_image.reduceRegion(
#             reducer=ee.Reducer.mean(),
#             geometry=geometry,
#             scale=30,
#             maxPixels=1e15,
#             bestEffort=True
#         ).getInfo()
        
#         annual_ndvi[year] = ndvi_value.get('NDVI', None)
    
#     return annual_ndvi

def plot_annual_ndvi_time_series(monthly_ndvi_results, output_path, start_year, end_year):
    """
    Plot monthly NDVI time series from monthly NDVI values for different categories and save to file.
    :param monthly_ndvi_results: List of dictionaries containing monthly NDVI values for each category.
    :param output_path: Path to save the plot image.
    :param start_year: Start year for the NDVI time series.
    :param end_year: End year for the NDVI time series.
    """
    plt.figure(figsize=(14, 7))
    
    for result in monthly_ndvi_results: 
        keterangan_label = result['keterangan']
        annual_ndvi = result['annual_ndvi']
        
        # Extract the keys (month-year) and values (NDVI)
        time_points = list(annual_ndvi.keys())
        ndvi_values = [annual_ndvi[time_point] for time_point in time_points]
        
        plt.plot(time_points, ndvi_values, marker='o', linestyle='-', label=f"{keterangan_label}")
    
    plt.xlabel('Time (Month-Year)')
    plt.ylabel('NDVI')
    plt.title(f'Monthly NDVI Time Series ({start_year} - {end_year})')
    plt.xticks(rotation=45)  # Rotate the x-axis labels for better readability
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# def plot_annual_ndvi_time_series(annual_ndvi_results, output_path, start_year, end_year):
#     """
#     Plot annual NDVI time series from annual NDVI values for different categories and save to file.
#     :param annual_ndvi_results: List of dictionaries containing annual NDVI values for each category.
#     :param output_path: Path to save the plot image.
#     :param start_year: Start year for the NDVI time series.
#     :param end_year: End year for the NDVI time series.
#     """
#     plt.figure(figsize=(12, 6))
    
#     for result in annual_ndvi_results: 
#         keterangan_label = result['keterangan']
#         annual_ndvi = result['annual_ndvi']
        
#         years = list(annual_ndvi.keys())
#         ndvi_values = [annual_ndvi[year] for year in years]
        
#         plt.plot(years, ndvi_values, marker='o', linestyle='-', label=f"{keterangan_label}")
    
#     plt.xlabel('Year')
#     plt.ylabel('NDVI')
#     plt.title(f'Annual NDVI Time Series ({start_year} - {end_year})')
#     plt.xticks(years, [str(year) for year in years], rotation=45)  # Format year numbers
#     plt.legend()
#     plt.grid(True)
#     plt.tight_layout()
#     plt.savefig(output_path)
#     plt.close()
  

def plot_geometries(map_results, output_path, default_crs='EPSG:3857'):
    """
    Plot the geometries on a map with unique colors for each label and save the map to a file.
    :param map_results: List of dictionaries containing'keterangan', and 'geometries'.
    :param output_path: Path to save the map image.
    :param default_crs: The default CRS to use if none is provided. Default is 'EPSG:3857'.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    # Define a color map
    color_map = plt.get_cmap('tab20')  # Use a colormap with enough distinct colors
    
    # Initialize bounding box coordinates
    minx, miny, maxx, maxy = np.inf, np.inf, -np.inf, -np.inf
    
    # Track handles for legend and label to color mapping
    handles = []
    labels = []
    label_color_map = {}

    # Plot each geometry with its corresponding color
    for i, result in enumerate(map_results):
        geometries = result['geometries']  # This might be a GeoDataFrame, GeoSeries, or individual geometry 
        keterangan_label = result['keterangan']
        label = f"{keterangan_label}"  # Simplified label
        
        # Convert to GeoSeries if necessary
        if not isinstance(geometries, gpd.GeoSeries):
            geometries = gpd.GeoSeries([geometries])

        # Ensure the geometries have a CRS; set to default if not
        if geometries.crs is None:
            geometries = geometries.set_crs(default_crs)
        
        # Ensure geometries are in Web Mercator projection for OSM overlay
        geometries = geometries.to_crs(epsg=3857)  # Web Mercator projection
        
        # Update bounding box coordinates
        bounds = geometries.total_bounds
        minx = min(minx, bounds[0])
        miny = min(miny, bounds[1])
        maxx = max(maxx, bounds[2])
        maxy = max(maxy, bounds[3])
        
        # Plot each geometry with the color corresponding to its label
        color = label_color_map.get(label, color_map(len(label_color_map) / 20))
        label_color_map[label] = color
        geometries.plot(ax=ax, edgecolor='k', facecolor=color, label=label)
        
        # Calculate centroids and add labels
        centroids = geometries.geometry.centroid
        for centroid in centroids:
            x, y = centroid.x, centroid.y
            ax.text(x, y, label, fontsize=8, ha='right', color=color)  # Simplified label
        
        # Add to legend handles and labels
        if label not in labels:
            handles.append(plt.Line2D([0], [0], marker='o', color='w', label=label, 
                                     markerfacecolor=color, markersize=10))
            labels.append(label)
    
    # Set axis limits to encompass all geometries with an added margin
    margin = 0.1  # 10% margin
    x_range = maxx - minx
    y_range = maxy - miny
    ax.set_xlim(minx - margin * x_range, maxx + margin * x_range)
    ax.set_ylim(miny - margin * y_range, maxy + margin * y_range)
    
    # Add OpenStreetMap base map
    ctx.add_basemap(ax, crs='EPSG:3857', source=ctx.providers.OpenStreetMap.Mapnik, zoom=10)
    
    ax.set_title('Geometries with Unique Colors for Each Label')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    # Limit the number of labels in the legend to avoid clutter
    ax.legend(handles=handles[:20], labels=labels[:20], loc='upper right', bbox_to_anchor=(1.15, 1))
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()



# Function to calculate total area of geometries
def calculate_total_area(geometries):
    """
    Calculate the total area of the geometries in hectares.
    :param geometries: GeoDataFrame with geometries.
    :return: Total area in hectares.
    """
    
    projected_geometries = geometries.to_crs('EPSG:3857')

    projected_geometries['area_ha'] = projected_geometries['geom'].area / 10000
    total_area_ha = projected_geometries['area_ha'].sum()
    return total_area_ha

# Function to calculate the variation in NDVI values
def calculate_variation(ndvi_results):
    """
    Calculate the variance of NDVI values for each category.
    :param ndvi_results: List of dictionaries containing NDVI values for each category.
    :return: List of dictionaries with variance information.
    """
    variation_results = []
    for result in ndvi_results: 
        keterangan_label = result['keterangan']
        monthly_ndvi = result['monthly_ndvi']
        ndvi_values = [monthly_ndvi[month] for month in range(1, 13) if monthly_ndvi[month] is not None]
        
        if len(ndvi_values) > 1:
            variance = pd.Series(ndvi_values).var()
        else:
            variance = None
        
        variation_results.append({  'keterangan': keterangan_label, 'variance': variance})
    
    return variation_results

# Function to save tabular information to a CSV file
def save_tabular_information(total_area_ha, variation_results, team_info, output_path):
    """
    Save the total area and variance information to a CSV file.
    :param total_area_ha: Total area in hectares.
    :param variation_results: List of dictionaries with variance information.
    :param team_info: The official team that has provided the public information.
    :param output_path: Path to save the CSV file.
    """
    # Identify the area with the highest variation
    highest_variation = max(variation_results, key=lambda x: x['variance']) if any(v['variance'] is not None for v in variation_results) else None 
    area_with_highest_variation = highest_variation['keterangan'] if highest_variation else 'N/A'
    variance_of_highest_variation = highest_variation['variance'] if highest_variation else 'N/A'
    
    # Make inferences on the variance
    variance_inference = "N/A"
    if variance_of_highest_variation != 'N/A':
        if variance_of_highest_variation > 0.5:
            variance_inference = "High variance observed, suggesting significant changes over time."
        elif variance_of_highest_variation > 0.2:
            variance_inference = "Moderate variance observed, indicating some level of change over time."
        else:
            variance_inference = "Low variance observed, implying stable conditions over time."

    # Prepare the tabular data
    tabular_data = {
        'Metric': [
            'Total Mangrove Area (hectares)', 
            'Official Team Providing Public Information',  
            'Area with Highest Variation', 
            'Variance of Highest Variation', 
            'Inference on Variance'
        ],
        'Value': [
            f"{total_area_ha:.2f}", 
            team_info, 
            area_with_highest_variation, 
            variance_of_highest_variation,
            variance_inference
        ]
    }
    
    # Save the data to a CSV file
    df = pd.DataFrame(tabular_data)
    df.to_csv(output_path, index=False)


# Main function to process NDVI extraction, plot time series, and generate report
def process_and_generate_report(credentials_path, service_account_email, db_params, year_param, area, output_dir='results'):
    """
    Process NDVI extraction, plot time series, and generate a report.
    :param credentials_path: Path to the service account credentials JSON file.
    :param service_account_email: Service account email for Earth Engine authentication.
    :param db_params: Dictionary containing database connection parameters.
    :param query: SQL query string to fetch geometries from the database.
    :param output_dir: Directory to save the results.
    :param year: Year for NDVI calculation.
    """
    # Create results directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Initialize Earth Engine
    initialize_earth_engine(credentials_path, service_account_email)

    # Create database engine
    engine = create_db_engine(db_params)
    query = f"""
            SELECT 1 id, keterangan, ST_Union(ST_Transform(geometry, 4326)) as geom 
            FROM staging.linked_data_view 
            WHERE keterangan ILIKE '%%{area}%%'
            GROUP BY keterangan
            """
    # Fetch geometries from the database
    geometries = fetch_geometries(query, engine)

    # Calculate total area
    total_area_ha = calculate_total_area(geometries)

    ndvi_results = []
    ndvi_yearly_results = []
    map_results = []
    
    for _, row in geometries.iterrows():
        geom = row['geom'] 
        keterangan_label = row['keterangan']
        
        gee_geom = convert_geom_to_gee(geom)
        monthly_ndvi = calculate_monthly_ndvi(gee_geom, year_param[1])
        yearly_ndvi = calculate_annual_ndvi(gee_geom, year_param[0], year_param[1] )
        
        ndvi_results.append({ 
            'keterangan': keterangan_label,
            'monthly_ndvi': monthly_ndvi
        })

        ndvi_yearly_results.append({ 
            'keterangan': keterangan_label,
            'annual_ndvi': yearly_ndvi
        })

        map_results.append({ 
            'keterangan': keterangan_label,
            'geometries': geom
        })
    
    # Plot NDVI time series
    plot_ndvi_time_series(ndvi_results, os.path.join(output_dir, f'{year_param[1]}_ndvi_time_series.png'), f'{year_param[1]}')
    plot_annual_ndvi_time_series(ndvi_yearly_results, os.path.join(output_dir, f'{year_param[0]}-{year_param[1]}_ndvi_time_series.png'), f'{year_param[0]}', f'{year_param[1]}')
    
    # Plot and save geometries on a map
    plot_geometries(map_results, os.path.join(output_dir, 'geometries_map.png'))

    # Calculate NDVI variation
    variation_results = calculate_variation(ndvi_results)
    
    # Save tabular information
    save_tabular_information(total_area_ha, variation_results, 'Dinas Pertanahan dan Penataan Ruang Kota Balikpapan', os.path.join(output_dir, 'summary_report.csv'))

    print("Processing complete. Results saved in", output_dir)

# Example usage
db_params = {
    'dbname': 'april_spatial_db',
    'user': 'april',
    'password': '2030april',
    'host': 'localhost',
    'port': '5434'
}


year_param = ['2018', '2023']

credentials_path = 'auth\gcp-learn-368515-cadc2c3f885c.json'
service_account_email = 'testing-earth-engine@gcp-learn-368515.iam.gserviceaccount.com'

process_and_generate_report(credentials_path, service_account_email, db_params, year_param, 'mangrove', output_dir='results' )