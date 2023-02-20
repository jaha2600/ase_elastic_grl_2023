#%% import modules
import numpy as np
import rasterio as rio
import configparser
import pandas as pd
from dem_functions import *
def main():
    # get inputs from config files
    config_file = 'configs_dem.ini'
    config = configparser.ConfigParser()
    config.read(config_file)
    # parse configuration file variables for functions
    geotiff = config.get('COREGISTER', 'file')
    dem_dir = config.get('COREGISTER', 'dem_directory_coreg')
    dem_ending = config.get('COREGISTER', 'dem_ending')
    pc_path = config.get('COREGISTER', 'pc_path')
    dem_dir_1 = config.get('CLOUD_CLIPPING', 'dem_directory_clouds')
    ref_dem_1 = config.get('CLOUD_CLIPPING', 'reference_dem_clouds')
    cloud_diff_thresh = config.getfloat('CLOUD_CLIPPING', 'cloud_diff_threshold')
    dem_dir_2 = config.get('TILE_DEMS', 'dem_directory_tile')
    grid_size = config.getfloat('TILE_DEMS', 'tile_grid_size')
    ref_dem_2 = config.get('TILE_DEMS', 'reference_dem_tile')
    pickl_file = config.get('TILE_DEMS', 'pickle_file_tiles')
    tied_ref_outdir = config.get('TILE_DEMS', 'ref_file_outdir_location')
    dem_dir_3 = config.get('CARST_DHDT', 'dem_directory_for_carst')
    ref_geom_root = config.get('CARST_DHDT', 'reference_geometry_root')
    point_cloud = config.get('CARST_DHDT', 'point_cloud_for_carst')
    carst_dem_dir = config.get('CARST_DHDT', 'carst_dem_directory')
    # run functions
    gtif2xyz(geotiff)
    coreg_dem_strips(dem_dir, dem_ending, pc_path)
    remove_clouds(dem_dir_1,ref_dem_1,cloud_diff_thresh)
    tile_dems(dem_dir_2,grid_size)
    tile_ref_dem(ref_dem_2, grid_size,pickl_file,tied_ref_outdir)
    make_dhdt_inputs(dem_dir_3,ref_geom_root,point_cloud)
    run_carst_dhdt(carst_dem_dir)
if __name__ == '__main__':
    main()
