#!/usr/local/anaconda3/bin/python
"""
Jasmine Hansen, 2023
Run median filter over image in commandline
"""

import sys
import rasterio
from scipy.signal import medfilt
import os
import argparse

#code to make it command line executable 
def getparser():
    parser = argparse.ArgumentParser(description="Run median filter on input GeoTiff")
    parser.add_argument('file', type=str, help='File')
    parser.add_argument('window_size', type=int, help='Filter Window Size - Must Be Odd No.')
    return parser

parser = getparser()
args = parser.parse_args()

filename = args.file
win_size = args.window_size

# function that will run a median filter over the input geotiff 
def main():
    path_list = [filename]
    output_path_list = []
    for pathname in path_list:
        fname, ext = os.path.splitext(pathname)
        med_name = '_med{}'.format(win_size)
        output_path_list.append(fname + med_name + ext)
    
    
    # This loop runs the filter over the dems
    for inputname,outputname in zip(path_list, output_path_list):
        print(inputname, outputname)
    
        with rasterio.open(inputname) as src:
            array = src.read()
            profile = src.profile
    
        # apply a X by X median filter to each band
        filtered = medfilt(array, (1,win_size,win_size)).astype('float32')
    
        # Write to tif, using the same profile as the source
        with rasterio.open(outputname, 'w', **profile) as dst:
            dst.write(filtered)

if __name__ == "__main__":
    main()

