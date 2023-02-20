# Sensitivity of Elastic Crustal Deformation in the Amundsen Sea Embayment

## DEFINE BEDROCK MASKS
`make_bedrock_masks.txt` contains a Google Earth Engine script to create custom defined reference point masks using NASA MEaSUREs (Rignot et al., 2017), and the dh/dt grid of Smith et al. (2020)

## DEM PROCESSING 
Code to process Digital Elevation Model (DEM) strips from the Reference Elevation Model of Antarctica (Howat et al., 2019) and produce high resolution surface change grids using Version 1.0 of the Cryosphere and Remote Sensing Toolkit (CARST).

There are four key stages:
1. Coregistration to reference point cloud
2. Removal of clouds using elevation threshold
3. Tiling of DEM strips to allow segmented processing
4. Run weighted linear regression using the Cryosphere And Remote Sensing Toolkit v1 (REF). For information on CARST usage please see: 

The main code is located in `dem_functions.py`and consists of the following functions:

1. `gtif2xyz` Convert a geotiff reference file into an xyz point cloud removing nodata values:
```
gtif2xyz(geotiff_name)
```
2. `coreg_dem_strips` Coregister DEM strips from REMA or ArcticDEM (REF) using the newly created reference point cloud:
```
coreg_dem_strips(dem_directory, dem_filename_ending, pointcloud_path)
```
3. `remove_clouds` Remove any remaining clouds from coregistered DEM strips by comparison to a reference DEM mosaic. The threshold value (max acceptable difference value) can be manually defined:
```
remove_clouds(dem_directory,reference_dem,difference_threshold)
```
4. `tile_dems` Tile DEM strips to allow processing in chunks for subsequent stages. User defines a tile size in km and indexed subdirectories are created with clipped DEMs. Two scripts, one for main DEM strips and one for reference geometry as required to run CARST dh/dt (see REF)
``` 
tile_dems(dem_directory,grid_size_km)
tile_ref_dem(ref_dem, grid_size,pickl_file,outdir)
```
5. `make_dhdt_inputs` Create correctly formatted input files for CARST. Outputs the required configuration file for each individual tile directory
```
make_dhdt_inputs(dem_directory,ref_geom_name,point_cloud)
```
6. `run_carst_dhdt` Run version 1.0 of CARST dh/dt algorithm (https://github.com/whyjz/CARST/releases/tag/v1.0-alpha)
    run_carst_dhdt(carst_dem_dir)

Instructions

Edit `configs_dem.ini` file to set correct locations for directories and required input files for each stage

Run `run_dem_processing.py` in the command line to process files:
```
python run_dem_processing.py
```
## FIGURES
Scripts used to create figures 2 and 3 for paper 
```
figure_2_create.py and figure_3_create.py
```

## References


CARST: 
Zheng, W., Durkin, W. J., Melkonian, A. K., & Pritchard, M. E. (2021, March 9). Cryosphere And Remote Sensing Toolkit (CARST) v2.0.0a1 (Version v2.0.0a1). Zenodo. http://doi.org/10.5281/zenodo.3475693

DEM Strips from REMA: 
Howat, I. M., Porter, C., Smith, B. E., Noh, M.-J., & Morin, P. (2019). The Reference Elevation Model of Antarctica. The Cryosphere, 13(2), 665–674. https://doi.org/10.5194/tc-13-665-2019

ICESat 2 dhdt:
Smith, B., Fricker, H. A., Gardner, A. S., Medley, B., Nilsson, J., Paolo, F. S., et al. (2020). Pervasive ice sheet mass loss reflects competing ocean and atmosphere processes. Science, 367(6496), 1239-1242. https://doi.org/10.1126/science.aaz5845

MEaSUREs
Rignot, E., J. Mouginot, & Scheuchl, B. (2017). MEaSUREs InSAR-Based Antarctica Ice Velocity Map (Version 2). [Dataset] NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/D7GK8F5J8M8R. 