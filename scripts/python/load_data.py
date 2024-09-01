import os
import logging
import pandas as pd
import geopandas as gpd
from sqlalchemy import create_engine, text

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Define the file names and relative paths
GPKG_FILE_NAME = 'lu.gpkg'
CSV_FILE_NAME = 'lu.csv'

def check_file_exists(file_path):
    """Check if the file exists at the given path."""
    if not os.path.isfile(file_path):
        logging.error(f"File not found: {file_path}")
        raise FileNotFoundError(f"File not found: {file_path}")
    logging.info(f"File found: {file_path}")

def create_schema_if_not_exists(engine, schema_name):
    """Create schema in PostgreSQL if it does not exist."""
    with engine.connect() as conn:
        try:
            conn.execute(text(f"CREATE SCHEMA IF NOT EXISTS {schema_name};"))
            logging.info(f"Schema '{schema_name}' is ready.")
        except Exception as e:
            logging.error(f"Error creating schema '{schema_name}': {e}")
            raise

def create_table_if_not_exists(engine, create_table_query):
    """Create table if it does not exist."""
    create_table_sql = f"""
    {create_table_query}
    """
    with engine.connect() as connection:
        connection.execute(text(create_table_sql))
        logging.info(f"Table created or already exists.")
        
def load_geodataframe(file_path):
    """Load GeoDataFrame from the given file path."""
    try:
        gdf = gpd.read_file(file_path)
        logging.info("GeoDataFrame loaded successfully.")
        return gdf
    except Exception as e:
        logging.error(f"Error loading GeoDataFrame: {e}")
        raise

def validate_crs(gdf):
    """Ensure CRS is defined and log its information."""
    if gdf.crs is None:
        logging.error("Coordinate Reference System (CRS) is not defined in the GeoDataFrame.")
        raise ValueError("Coordinate Reference System (CRS) is not defined in the GeoDataFrame.")
    logging.info(f"Geometry column name: {gdf.geometry.name}")
    logging.info(f"Geometry column type: {gdf.geometry.dtype}")


def load_csv_to_dataframe(file_path):
    """Load CSV file into a DataFrame."""
    try:
        df = pd.read_csv(file_path)
        logging.info("CSV file loaded into DataFrame successfully.")
        return df
    except Exception as e:
        logging.error(f"Error loading CSV file: {e}")
        raise

def add_id_column(df, id_column_name='id'):
    """Add a sequential ID column to the DataFrame or GeoDataFrame and reorder columns to place the ID column first."""
    # Add ID column
    df[id_column_name] = range(1, len(df) + 1)
    
    # Reorder columns to place ID column first
    columns = [id_column_name] + [col for col in df.columns if col != id_column_name]
    df = df[columns]
    
    return df


def upload_geodataframe_to_postgis(gdf, table_name, engine, schema):
    """Upload GeoDataFrame to PostgreSQL with the given table name and schema, adding an ID column."""
    try:
        # Add ID column to GeoDataFrame
        gdf = add_id_column(gdf)
        
        # Upload GeoDataFrame to PostgreSQL
        gdf.to_postgis(name=table_name, con=engine, schema=schema, if_exists='replace', index=False)
        logging.info(f"GeoDataFrame loaded successfully into PostgreSQL table '{table_name}'.")
    except Exception as e:
        logging.error(f"Error uploading GeoDataFrame to PostgreSQL: {e}")
        raise

def upload_dataframe_to_postgres(df, table_name, engine, schema):
    """Upload DataFrame to PostgreSQL with the given table name, adding an ID column."""
    try:
        # Add ID column to DataFrame
        df = add_id_column(df)
        
        # Upload DataFrame to PostgreSQL
        df.to_sql(name=table_name, con=engine, schema=schema, if_exists='replace', index=False)
        logging.info(f"DataFrame loaded successfully into PostgreSQL table '{table_name}'.")
    except Exception as e:
        logging.error(f"Error uploading DataFrame to PostgreSQL: {e}")
        raise

def main():
    # Define the relative paths to the data files
    script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    geopackage_file_path = os.path.join(script_dir, '..', 'data', GPKG_FILE_NAME)  # GeoPackage file path
    csv_file_path = os.path.join(script_dir, '..', 'data', CSV_FILE_NAME)  # CSV file path

    # Define database connection and schema
    engine = create_engine('postgresql://april:2030april@localhost:5434/april_spatial_db')
    schema = 'staging'
    geotable_name = 'tb_lu_dataset'
    csv_table_name = 'tb_lu_csv_dataset' 

    try:
        # Check if the GeoPackage file exists
        check_file_exists(geopackage_file_path)

        # Create schema if not exists
        create_schema_if_not_exists(engine, schema)

        # Load GeoDataFrame
        gdf = load_geodataframe(geopackage_file_path)

        # Validate CRS
        validate_crs(gdf)

        # Upload GeoDataFrame to PostgreSQL
        upload_geodataframe_to_postgis(gdf, geotable_name, engine, schema)

        # Check if the CSV file exists
        check_file_exists(csv_file_path)

        # Load CSV into DataFrame
        df = load_csv_to_dataframe(csv_file_path)
 
        # Upload DataFrame to PostgreSQL
        upload_dataframe_to_postgres(df, csv_table_name, engine, schema)

    except Exception as e:
        logging.error(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
