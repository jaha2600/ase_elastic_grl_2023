import subprocess
from osgeo import gdal
import os, shutil, sys
import numpy as np
from pygeotools.lib import iolib, warplib
from osgeo import gdal
import os, shutil, math, sys
import rasterio
import numpy as np
from rasterio.mask import mask
from shapely.geometry import box
import geopandas as gpd
from fiona.crs import from_epsg
import pycrs
import pickle
import warnings
import pandas as pd

def writeFile(filename,rasterOrigin,pixelWidth,pixelHeight,geoprojection,data):
    cols = data.shape[1]
    rows = data.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]
    
    #(x,y) = data.shape
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    # you can change the dataformat but be sure to be able to store negative values including -9999
    dst_datatype = gdal.GDT_Float32
    dst_ds = driver.Create(filename,cols,rows,1,dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    dst_ds.SetProjection(geoprojection)
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)

def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

def gtif2xyz(geotiff):
    # convert to xyz file 
    xyz_outfile = os.path.splitext(geotiff)[0] + '.xyz'
    command='gdal2xyz.py -band 1 {} > {}'.format(geotiff, xyz_outfile)
    subprocess.run(command,shell=True)
    # remove lines with nodata val
    sed_cmd="sed '/nan/d' -i {}".format(xyz_outfile)

def coreg_dem_strips(dem_dir, dem_ending, pc_path):
    #paths for ames stereo pipeline and dem coregistration
    ASP_CODE='/usr/local/StereoPipeline/bin/'
    CODE='/home/jasmine/Applications/demcoreg/demcoreg/apply_dem_inv_translation.py'
    CSVDATA = pc_path
    PC_NAME=os.path.split((os.path.splitext(CSVDATA)[0])[1])
    # state the projection of the point cloud (string below is for EPSG:3031)
    CSVPROJECTION="'+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'"
    #set the format of the point cloud file (x,y,z)
    CSVFORMAT="'1:easting 2:northing 3:height_above_datum'"
    DEM_PATH = dem_dir
    path_list = [os.path.join(DEM_PATH,fname) for fname in os.listdir(DEM_PATH) if fname.endswith(dem_ending)]
    os.chdir(DEM_PATH)

    #for each file in list run pc_align 
    #produces a new subdirectory called CORRECTED_point_cloud_name, within this are the output files from the pc_align algorithm.
    #every pc_align run produces a *pc_align*.txt file, regardless of whether or not it was successful or not.
    #*trans_reference.tif files are ONLY produced if the pc_align run is successful
    print('Running Step 1 - PC Align')
    for f in path_list:
        print('Running PC_ALIGN on {}'.format(f))
        name=os.path.split(os.path.splitext(f)[0])[1]
        total_name = name + '.tif'
        ASP_CMD ="{}pc_align --max-displacement 50 --tif-compress=NONE --save-inv-transformed-reference-points --threads 32 -o CORRECTED/{} --csv-proj4 '{}' --csv-format '{}' {} $CSVDATA".format(ASP_CODE,CSVPROJECTION,CSVFORMAT,name,total_name,CSVDATA)
        subprocess.run(ASP_CMD,shell=True)
    print('Stage 1 Complete')

    COR_PATH=os.path.join(DEM_PATH,'CORRECTED')
    os.chdir(COR_PATH)
    # list the trans reference.tif files which will show which files the pc_align algorithm was successful on
    # save only the root of the file (i.e. file id, date, satellite, segment number etc.)
    list_transl='ls *trans_reference.tif | cut -d"." -f1 > trans_list'
    subprocess.run(list_transl,shell=True)
    root_list='cat trans_list | cut -d"_" -f1-5 > trans_root_list'
    subprocess.run(root_list,shell=True)
    #check there is no existing pc_align list and if there is remove it, otherwise it just appends the new ones to the bottom
    pc_filelist_name = 'pc_file_list_{}'.format(os.path.split(os.path.splitext(CSVDATA)[0])[1])
    if os.path.exists(pc_filelist_name) == True:
        rm_file = 'rm {}'.format(pc_filelist_name)
        subprocess.run(rm_file,shell=True)
    else:
        print('Processing')
    ##########
    print('Identifying succesful pc_align runs')
    #list pc_align files for successful pc_align runs
    cmd1 = "for file_root in $(cat trans_root_list) ; do ls $PWD/${file_root}*pc_align-*.txt >> pc_file_list_{} ; done".format(PC_NAME)
    # CHECK THIS LINE
    cmd2 = 'for file_root in $(cat trans_root_list) ; do cp pc_file_list_{} {} ; done'.format(PC_NAME,DEM_PATH)
    subprocess.run(cmd1,shell=True)
    subprocess.run(cmd2,shell=True)

    #move to main directory
    os.chdir(DEM_PATH)
    print('Running Step 2 - Apply Inverse Transform to DEM')
    # for each pc align file run apply_dem_translation.py 
    #this takes the inverse transform from the pc align file and applies it
    # OPEN pc filelist and use that in the loop 
    pc_filelist = open(pc_filelist_name, 'r')
    pc_file_data = pc_filelist.read()
    pc_into_list = pc_file_data.split("\n")
    pc_filelist.close()
    for pc in pc_into_list:
        print('Applying Translation to {}'.format(pc))
        dem_root = pc.split('-')[0]
        dem_filename = os.path.join(dem_root,'tif')
        dem_filename_shean='{}_{}_trans.tif'.format(dem_root,PC_NAME)
        # apply inverse translation on each DEM file
        trans_cmd='python {} {} {} {}'.format(CODE,dem_filename,pc,PC_NAME)
        subprocess.run(trans_cmd,shell=True)
        # resample DEM to 30m 
        gdal_cmd = 'gdal_warp -tr 30 30 -r bilinear {} {}_{}_trans_30.tif'.format(dem_filename_shean,dem_root,PC_NAME)
        subprocess.run(gdal_cmd,shell=True)   
        # run compression over the 30 m files        
        m_com="gdal_translate -co 'COMPRESS=LZW' -co bigtiff=if_safer {}_{}_trans_30.tif {}_{}_trans_30m.tif".format(dem_root,PC_NAME,dem_root,PC_NAME)
        subprocess.run(m_com,shell=True)
        o_com="gdal_translate -co 'COMPRESS=LZW' -co bigtiff=if_safer {}_{}_trans.tif {}_{}_trans.tif".format(dem_root,PC_NAME,dem_root,PC_NAME)
        subprocess.run(o_com,shell=True)
        # run compression over the native res files
    
    print('Step 2 Complete')
   
def remove_clouds(dem_dir_1,ref_dem,thresh):
   # parser.add_argument('dem_dir', type=str, help='DEM Strip directory')
   # parser.add_argument('ref_dem', type=str, help='Reference Dem')
   # parser.add_argument('thresh', type=int, help='Difference Threshold (i.e. 10 = -10 to 10')
    DEM_DIRECTORY = dem_dir_1
    ref_tile = ref_dem
    threshold = thresh
    if not iolib.fn_check(DEM_DIRECTORY): 
        sys.exit("Unable to find DEM strip directory: %s" % DEM_DIRECTORY)
    if not iolib.fn_check(ref_tile):
        sys.exit("Unable to find Reference DEM: %s" % ref_tile)
    #remove clouds by comparing to a reference dem. 
    coreg_dir_path = DEM_DIRECTORY
    cloud_clip_path_list = [os.path.join(coreg_dir_path, fname) for fname in os.listdir(coreg_dir_path) if fname.endswith("trans_30m.tif")]
    
    #threshold values for masking - everything between these values is retained.
    min_difference_value = (-(threshold))
    max_difference_value = (threshold)
    
    #function to convert numpy array to geotiff.
    for i, dem in enumerate(cloud_clip_path_list):
        print(i, "/", (len(cloud_clip_path_list) -1))
        #use pygeotools to make rema the same extent
        dem_fn_list = [ref_tile, dem]
        ds_list = warplib.memwarp_multi_fn(dem_fn_list, extent='intersection', res='min', t_srs=dem)
    
        #open raw dem in gdal to get specifications.
        dataset = gdal.Open(dem)
        #information about the raster
        cols = dataset.RasterXSize 
        rows = dataset.RasterYSize
        projection = dataset.GetProjection()
        geotransform = dataset.GetGeoTransform()
        xMin = geotransform[0]
        yMax = geotransform[3]
        xMax = xMin + (dataset.RasterXSize/geotransform[1])
        yMin = yMax + (dataset.RasterXSize/geotransform[5])
        dataset_extent = (xMin, xMax, yMin, yMax)
        print('extent: ',dataset_extent)
        #get data
        data_raster = dataset.GetRasterBand(1)
        #get no data value
        nodata_value = data_raster.GetNoDataValue()
        #Convert raster to array
        dataset_array = dataset.GetRasterBand(1).ReadAsArray(0,0,cols,rows).astype(np.float)
        #Load datasets to NumPy masked arrays, part of this is redundant need to change
        rema_dem, input_dem_file  = [iolib.ds_getma(i) for i in ds_list]
        #Calculate elevation difference between rema and dem strip
        test_diff = rema_dem - input_dem_file
        #set if statements to generate a mask for the image that is 1 for good data and 0 for bad
        test_diff[test_diff < min_difference_value] = 0
        test_diff[test_diff > max_difference_value] = 0
        test_diff[test_diff != 0 ] = 1
        # generate the new filename
        output_section = os.path.split(dem)[1]
        # multiply input dem with the new mask to get new dem
        clipped_dem = input_dem_file * test_diff
        #set zeros to nodata value 
        clipped_dem[clipped_dem == 0.0] = nodata_value
        #save out this file
        output_filename_new = os.path.join(os.path.splitext(output_section)[0] + '_clip.tif') 
        new_name_final = os.path.join(coreg_dir_path, output_filename_new)
        # save out new filename
        writeFile(new_name_final,(xMin,yMax),30,-30,projection,np.array(clipped_dem,dtype=float))
    if not os.path.exists(os.path.join(coreg_dir_path + 'cloud_masked')):
        os.mkdir(os.path.join(coreg_dir_path + 'cloud_masked'))

    clipped_pathlist = [os.path.join(coreg_dir_path, fname) for fname in os.listdir(coreg_dir_path) if fname.endswith("trans_30m_clip.tif")]
    #move cloud masked files into their own directory         
    for a in clipped_pathlist:
        init_loc = coreg_dir_path + (os.path.basename(a)) 
        end_loc = coreg_dir_path + "cloud_masked/" + (os.path.basename(a))
        shutil.move(init_loc, end_loc)

def tile_dems(dem_dir_2,grid_size):
    warnings.simplefilter(action='ignore', category=FutureWarning)
    #parser.add_argument('dem_dir', type=str, help='DEM Strip directory')
    #parser.add_argument('grid_size', type=int, help='Grid Size in Km')
    DEM_DIRECTORY = dem_dir_2
    grid_size_km = grid_size
    if not iolib.fn_check(DEM_DIRECTORY): 
        sys.exit("Unable to find DEM strip directory: %s" % DEM_DIRECTORY)

    coreg_dir_path = DEM_DIRECTORY
    directory = coreg_dir_path
    grid_size = grid_size_km
    GRID_BOX_RESOLUTION_M = grid_size * 1000

    path_list = [os.path.join(directory, fname) for fname in os.listdir(directory) if fname.endswith("_clip.tif")]
    
    # List of bounding box coordinates for each file (above)
    x_min_list = [None] * len(path_list)
    x_max_list = [None] * len(path_list)
    y_min_list = [None] * len(path_list)
    y_max_list = [None] * len(path_list)
    
    x_resolution = None
    y_resolution = None
    
    #names for the directoriees
    x_range = list(range(0,11))
    y_range = list(range(0,9)) 
            
    ## 1. Loop through all the files, get a list of bounding-box coordinates.
    for i,path in enumerate(path_list):
        info = gdal.Info(path, format="json")
    #    print(info)
        corner_coords = info['cornerCoordinates']
        upper_left_coords = corner_coords['upperLeft']
        upper_right_coords = corner_coords['upperRight']
        lower_left_coords = corner_coords['lowerLeft']
        lower_right_coords = corner_coords['lowerRight']
    
        file_x_res, file_y_res = info['geoTransform'][1], info['geoTransform'][5]
        # Get the pixel resolution from the first file we read.
        if x_resolution is None and y_resolution is None:
            x_resolution = file_x_res
            y_resolution = file_y_res
        # Assert that the resolution of this file is identical to the previous files
        else:
            assert((x_resolution == file_x_res) and (y_resolution == file_y_res))
        
    #    print(upper_left_coords, upper_right_coords, lower_left_coords, lower_right_coords)
    
        # Fill in the x_min value
        x_min_list[i] = min(upper_left_coords[0], lower_left_coords[0])
        x_max_list[i] = max(upper_right_coords[0], lower_right_coords[0])
        y_min_list[i] = min(lower_left_coords[1], lower_right_coords[1])
        y_max_list[i] = max(upper_left_coords[1], upper_right_coords[1])
    
    
    grid_x_min = min(x_min_list)
    grid_y_min = min(y_min_list)
    grid_x_max = max(x_max_list)
    grid_y_max = max(y_max_list)
    
    print("Final grid measures {0} by {1} km.".format((grid_x_max - grid_x_min)/1000.,
                                                      (grid_y_max - grid_y_min)/1000.))
    
    grid_size_x = (grid_x_max - grid_x_min)
    grid_size_y = (grid_y_max - grid_y_min)
    
    #work out how many subset boxes you need rounded up to nearest whole number
    x_cell_number = math.ceil(grid_size_x/GRID_BOX_RESOLUTION_M)
    y_cell_number = math.ceil(grid_size_y/GRID_BOX_RESOLUTION_M)
    
    print("Subset grid will be {0} boxes across and {1} boxes down".format((x_cell_number),(y_cell_number)))
    
    # make a list of x values for the subset grids.
    grid_x_min_list = [0] * x_cell_number
    #assign first value
    grid_x_min_list[0] = grid_x_min
    #do the same with y values
    grid_y_min_list = [0] * y_cell_number
    grid_y_min_list[0] = grid_y_min
    

    
    #loop through the x_list and generate the minimum x coordinates for each subset box
    #these points are the 'bottom left'
    for index in range(1,len(grid_x_min_list)):
        grid_x_min_list[index] = grid_x_min_list[index-1]+GRID_BOX_RESOLUTION_M
    
    #do the same for the y coordinate.
    for ind in range(1,len(grid_y_min_list)):
        grid_y_min_list[ind] = grid_y_min_list[ind-1] + GRID_BOX_RESOLUTION_M
     
    #fname = 'berp_outline_new.tif'
    #new_directory = '/data2/ANTARCTICA_2019/SHAPEFILES/POLENET_area_outlines_for_dhdt/'
    #new_path_list = [os.path.join(new_directory, fname)]
    
    #for the first dem in the list find what indices it lies between in thte x and y arrays.
    
    for i, grid_x_min_value in enumerate(grid_x_min_list):
        for j, grid_y_min_value in enumerate(grid_y_min_list):
            print('Working on grid cell {}_{}'.format(i,j))
            for img_i, dem_name in enumerate(path_list):
                # If no part of the image bounding box overlaps the grid box, then
                # skip this image and move along.
                if not ((x_min_list[img_i] <= grid_x_min_value + GRID_BOX_RESOLUTION_M) and \
                        (x_max_list[img_i] >= grid_x_min_value) and \
                        (y_min_list[img_i] <= grid_y_min_value + GRID_BOX_RESOLUTION_M) and \
                        (y_max_list[img_i] >= grid_y_min_value)):
                    continue
                #print(i)
                #print(j)
                
                file = rasterio.open(dem_name)    
                #out_filename = os.path.join(os.path.splitext(path_list[0]))
                
                minx, miny = grid_x_min_value, grid_y_min_value
                maxx, maxy = (grid_x_min_value + GRID_BOX_RESOLUTION_M), (grid_y_min_value + GRID_BOX_RESOLUTION_M)
                bounding_box = box(minx, miny, maxx, maxy)
                
                ####
                geo = gpd.GeoDataFrame({'geometry': bounding_box}, index = [0], crs=from_epsg(3031))
                geo = geo.to_crs(crs=file.crs.data['init'])
        
                coords = getFeatures(geo)
        
                out_img, out_transform = mask(file, coords, crop=True)
                
                out_meta = file.meta.copy()
                
                #parse EPSG value from CRS so can make proj4 string
                #epsg_code = int(file.crs.data[5:])
                epsg_code = int(file.crs.data['init'][5:])
                #print(epsg_code)
        
    #update the metadata with dimensions, transform and crs
                out_meta.update({"driver": "GTiff",
                                 "height": out_img.shape[1],
                                 "width": out_img.shape[2],
                                 "transform": out_transform,
                                 "crs": pycrs.parse.from_epsg_code(epsg_code).to_proj4(),
                                 "nodata": -9999,
                                 })
        
                #generate the new filename
                
                #original dem filenam
                #new_dem_name = [os.path.splitext(dem_name)[0]]
                #path_name_only
                file_string ="{}_{}_{}_{}.tif"
               
                new_file_path = (file_string.format(dem_name[0:len(dem_name)-4], i, j, grid_size))
                # this doesnt apppear to be working in the correct way
                # convert out image max value to integer for if statement
                #if out_img.max() != "nan": 
                    #out_image_max_int = int(out_img.max())
                #if the out image max value is very low then the dem only contains nan (-9999) data so do not save it out.
                
                if np.isnan(out_img.max()) == True:
                    continue
                elif ((int(out_img.max()) < (-4000))):
                    continue
                else:
                    with rasterio.open(os.path.join(new_file_path), "w", **out_meta) as dest:
                        dest.write(out_img)

    # save grid_x_min and grid_y_min lists to a pickle file so they can be used when tiling reference_dems. 
    
    pickle_name = directory + 'grid_lists.pkl'
    with open(pickle_name, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([grid_x_min_list, grid_y_min_list], f)

    
    for dir_1, dir_1_value in enumerate(grid_x_min_list):
        for dir_2, dir_2_value in enumerate (grid_y_min_list):
            for path_1, path_1_value in enumerate (path_list):
                
                directory_name = "{}{}_{}_{}/"
                new_directory_name = (directory_name.format(directory, dir_1, dir_2,grid_size))
                try:
                    os.mkdir(new_directory_name)
                except OSError:
                    pass
                #os.mkdir(new_directory_name)
                
                #make new path to file
                split_path_names = os.path.split(path_1_value)
                split_file_name = (os.path.splitext(split_path_names[1])[0])
                
                current_new_dem_structure = "{}{}_{}_{}_{}.tif"
                current_location_dem_subset_filename = (current_new_dem_structure.format(directory, split_file_name, dir_1, dir_2,grid_size))
                dem_subset_filename = os.path.split(current_location_dem_subset_filename)[1]
                
                grid_id_part = "{}_{}_{}"
                grid_id = grid_id_part.format(dir_1,dir_2, grid_size)
                
                
                new_location_subset_filename = (os.path.join(directory, grid_id + "/", dem_subset_filename))
                
                #if the file exists then move it into a new directory
                if os.path.exists(current_location_dem_subset_filename) == True:
                    shutil.move(current_location_dem_subset_filename, new_location_subset_filename)
                else:
                    continue
    #remove empty directories at the end.
    folders = list(os.walk(directory))[1:]
    for folder in folders:
        if not folder[2]:
            os.rmdir(folder[0])

def tile_ref_dem(ref_dem, grid_size,pickl_file,outdir):
#for the geopandas variable geo Future Warning
    warnings.simplefilter(action='ignore', category=FutureWarning)
    #function to tile coregistered, cloud clipped dems into tiles.
    #parser.add_argument('ref_dem', type=str, help='Reference DEM to Tile')
    #parser.add_argument('grid_size', type=int, help='Grid Size in Km')
    #parser.add_argument('pickl_file', type=str, help='Pickle File for Run to get extents')
    #parser.add_argument('-outdir', default=None, help='Output directory')
   

    dhdt_ref_dem = ref_dem
    grid_size_km = grid_size
    pickle_file = pickl_file
    outdir = outdir

    path_root = os.path.split(dhdt_ref_dem)[0]
    outfile_path = path_root + '/' + outdir

    if not iolib.fn_check(dhdt_ref_dem): 
        sys.exit("Unable to find ref DEM: %s" % dhdt_ref_dem)
    if not iolib.fn_check(pickle_file): 
        sys.exit("Unable to find ref DEM: %s" % pickle_file)
    if not iolib.fn_check(outfile_path):
        os.makedirs(outfile_path)
    
    directory = pickle_file
    dem_name = dhdt_ref_dem
    
    with open(directory, 'rb') as f:  # Python 3: open(..., 'rb')
        grid_x_min_list, grid_y_min_list = pickle.load(f) 
    
    grid_size = grid_size_km
    GRID_BOX_RESOLUTION_M = grid_size * 1000 
    
    for i, grid_x_min_value in enumerate(grid_x_min_list):
        for j, grid_y_min_value in enumerate(grid_y_min_list):
                print('Working on grid cell {}_{}'.format(i,j))
                file = rasterio.open(dem_name)    
                #out_filename = os.path.join(os.path.splitext(path_list[0]))
                
                minx, miny = grid_x_min_value, grid_y_min_value
                maxx, maxy = (grid_x_min_value + GRID_BOX_RESOLUTION_M), (grid_y_min_value + GRID_BOX_RESOLUTION_M)
                bounding_box = box(minx, miny, maxx, maxy)
                
                ####
                geo = gpd.GeoDataFrame({'geometry': bounding_box}, index = [0], crs=from_epsg(3031))
                geo = geo.to_crs(crs=file.crs.data)
        
                coords = getFeatures(geo)
        
                out_img, out_transform = mask(file, coords, crop=True)
                
                out_meta = file.meta.copy()
                
                #parse EPSG value from CRS so can make proj4 string
                epsg_code = int(file.crs.data['init'][5:])
                #print(epsg_code)
        
                #update the metadata with dimensions, transform and crs
                out_meta.update({"driver": "GTiff",
                                 "height": out_img.shape[1],
                                 "width": out_img.shape[2],
                                 "transform": out_transform,
                                 "crs": pycrs.parse.from_epsg_code(epsg_code).to_proj4(),
                                 "nodata": -9999,
                                 })
        
                #generate the new filename
                
                #path_name_only
                file_string ="{}{}_{}_{}_{}.tif"
                new_filename = '/ref_geometry_' + outdir
                new_file_path = (file_string.format(outfile_path, new_filename, i, j, grid_size))

                if np.isnan(out_img.max()) == True:
                    continue
                elif ((int(out_img.max()) < (-4000))):
                    continue
                else:
                    with rasterio.open(os.path.join(new_file_path), "w", **out_meta) as dest:
                        dest.write(out_img)           


def make_dhdt_inputs(dem_dir,ref_geom_root,point_cloud):
    warnings.simplefilter(action='ignore', category=FutureWarning)

#example inputs:
#direct = '/data/ANTARCTICA/ORIGINAL_DEMS/all_dems_sep_2020/8m/trans_30m/good/'
#site_id = 'ase_loose'
#ref_root = ref_root = '/data/ANTARCTICA/files_for_dhdt/ase_sep_2020/reference_geometry_files/geotiffs/dotson_ref_geometry'
# get information
   # parser.add_argument('dem_dir', type=str, help='DEM Directory')
   # parser.add_argument('ref_geom_root', type=str, help='Path and filename for ref geometry')
   # parser.add_argument('point_cloud', type=str, help='pointcloud name')

    DEM_DIRECTORY = dem_dir
    ref_geotiff_root = ref_geom_root
    site_id = point_cloud

    if not iolib.fn_check(DEM_DIRECTORY): 
        sys.exit("Unable to find ref DEM: %s" % DEM_DIRECTORY)

# generate necessary input files for the dhdt script 
    coreg_dir_path = DEM_DIRECTORY
    directory =coreg_dir_path
    
    
    #import the masterfile into python as a pandas dataframe, add a date column, save back out as csv?  
    uncertainties_master_file = pd.read_csv(coreg_dir_path + (("{}_master_uncertainties_30m.csv").format(site_id)), usecols=[0,1,2,3])
    header_name = ["Filename", "Date", "Uncertainty", "Offset to Ref Points"]
    uncertainties_master_file.columns=header_name
    new = uncertainties_master_file["Filename"].str.split("/", expand = True)
    uncertainties_master_file["dates"] = new.iloc[:,-1].str[:8]
    uncertainties_master_file["dates"] = pd.to_datetime(uncertainties_master_file.dates)
    uncertainties_master_file['Date'] = uncertainties_master_file['dates'].dt.strftime("%Y-%m-%d")
    del uncertainties_master_file["dates"]
    #make a column for the paths alone - change as needs to be all but end.......
    # drop last column in new so its only the path
    uncertainties_master_file['Master_Filename'] = new.iloc[:,-1]
    new.drop(new.columns[len(new.columns)-1], axis=1, inplace=True)
    #concatenate all of those columns together
    uncertainties_master_file['Master_Paths'] = new.apply(lambda x: "/".join(x.dropna().astype(str).values), axis=1)
    uncertainties_master_file['Master_Paths'] = uncertainties_master_file['Master_Paths'].astype(str) + "/"
    
    uncertainties_master_file['Master_Filename'] = uncertainties_master_file['Master_Filename'].str.replace(r'.tif$', '')
    #ventually put this in a loop. 
    
    tile_filelist = [os.path.join(directory, fname) for fname in os.listdir(directory) if fname.endswith("_70")]
    
        
    for tile_directory in tile_filelist:
        print(tile_directory)
        tile_dem_list = [os.path.join(directory, file) for file in os.listdir(tile_directory) if file.endswith('70.tif')]
        #test_list = []
        uncertainties_new_file = pd.DataFrame()
        for item in tile_dem_list:
            value = item.split("/")[-1]
            value_2 = value.split('_')[:-4]
            joined = '_'.join(value_2)
            #test_list.append(joined)
            uncertainties_new_file = uncertainties_new_file.append(uncertainties_master_file[uncertainties_master_file.Master_Filename.str.contains(joined,case=False)])
        uncertainties_new_file['NewPathFile'] = tile_directory + "/" + uncertainties_new_file['Master_Filename'] + "_clip_" + (os.path.split(tile_directory)[1]) + ".tif"
        del uncertainties_new_file['Filename']
        del uncertainties_new_file['Master_Paths']
        del uncertainties_new_file['Master_Filename']
        cols = uncertainties_new_file.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        uncertainties_new_file = uncertainties_new_file[cols]
        original_headers = ['filename', 'date', 'uncertainty', 'mean_offset_wrt_refpts']
        uncertainties_new_file.columns=original_headers
        #save out the file. 
        uncertainty_filename = tile_directory + "/" + ((os.path.split(tile_directory)[1]) + "_" + "uncertainties.txt")
        uncertainties_new_file.to_csv(uncertainty_filename, index = False)
    
                                       
        # now make the default file for each tile directory. 
        #need to have the correct filenames and structure.
        reference_geotiff = (ref_geotiff_root + '_{}.tif').format((os.path.split(tile_directory)[1]))
        
        file = open(tile_directory + '/' + (os.path.split(tile_directory)[1]) + "_default.ini",'w')
        file.write('[demlist]\n')
        file.write('# === DEM List ===\n')
        csv_file = ('csvfile = {}/{}_uncertainties.txt').format(tile_directory, os.path.split(tile_directory)[1])
        file.write(csv_file)
        file.write('\n')
        file.write('\n')
        file.write('[refgeometry]\n')
        file.write('# ==== Reference Geometry, as a GeoTiff file ====\n')
        file.write('# ==== The grid spacing and the map extent will be used. ====\n')
        file.write('# ==== Note that the SRS must be same with all the input DEMs! ====\n')
        gtiff_name = ('gtiff = {}').format(reference_geotiff)
        file.write(gtiff_name)
        file.write('\n')
        file.write('\n')
        file.write('[settings]\n')
        refdate_name = (('refdate = {}').format(uncertainties_new_file['date'].min()))
        file.write(refdate_name)
        file.write('\n')
        file.write('\n')
        file.write('[result]\n')
        file.write('# ==== DHDT Result Options ====\n')
        picklefile = ('picklefile = {}/refgeo_30m_{}_{}.p').format(tile_directory, site_id, (os.path.split(tile_directory)[1]))
        file.write(picklefile)
        file.write('\n')
        dhdt_name = ('dhdt_prefix = {}/{}_30m_{}').format(tile_directory, site_id, (os.path.split(tile_directory)[1]))
        file.write(dhdt_name)
        file.close()

def run_carst_dhdt(dem_dir_3):
    DEM_DIRECTORY = dem_dir_3
    coreg_dir_path = DEM_DIRECTORY
    directory = coreg_dir_path 
    tile_directory_filelist = [os.path.join(directory, fname) for fname in os.listdir(directory) if fname.endswith("_70")]
    for i, dem_tiles in enumerate(tile_directory_filelist):
        print(i+1, "/", (len(tile_directory_filelist)))
        print('Running dhdt code for ' + (os.path.split(dem_tiles)[1]))
        dhdt_command = ('python /home/jasmine/code/CARST-master/dhdt/dhdt.py ' + dem_tiles + "/" + (os.path.split(dem_tiles)[1]) +"_default.ini")
        os.system(dhdt_command)