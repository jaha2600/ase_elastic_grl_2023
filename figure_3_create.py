#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Jasmine Hansen, 2023
"""

import matplotlib.pyplot as plt
import pickle
import numpy as np
import pandas as pd
import os
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# Then igfs is a dictionary keyed by the H3 resolution. 
# so igfs[6], igfs[7], igfs[8], and igfs[9]. 

# These all have keys 'dist_deg', 'dist', and 'dUdH'
# And the ones you will need are dUdH on the y and dist (in km) on the x.

# Let me know if you have any questions

file = 'igf_disk_experiment.pkl'
#filelist = [f for f in os.listdir(directory) if f.endswith('.pkl')]



# filename = file
# parts = filename.split('_')
# file_ending = parts[len(parts)-1].split('.')[0]
   
igfs = pickle.load(open(file,'rb'))

r6_dist = igfs[6]['dist']
r7_dist = igfs[7]['dist']
r8_dist = igfs[8]['dist']
r9_dist = igfs[9]['dist']

r6_dudh = igfs[6]['dUdH']
r7_dudh = igfs[7]['dUdH']
r8_dudh = igfs[8]['dUdH']
r9_dudh = igfs[9]['dUdH']



#creating a array of values between
#0 to 10 with a difference of 0.1
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
#creating figure and axes object
fig, ax = plt.subplots()

#plotting the curve
ax.plot(r6_dist, r6_dudh,c = CB_color_cycle[0],lw=2)
ax.plot(r7_dist, r7_dudh, c = CB_color_cycle[1],lw=2)
ax.plot(r8_dist, r8_dudh,c = CB_color_cycle[2],lw=2)
ax.plot(r9_dist, r9_dudh,c = CB_color_cycle[5],lw=2)
#formatting axes
ax.set_xlabel("Distance from Origin Point (km)")
ax.set_ylabel("Vertical Displacement (mm)")
#ax.set_title("Sine Wave") 
leg = '0.175 km', '0.46 km', '1.22 km', '3.23 km'
ax.legend(leg, title='Disc Radius', loc='lower right')
ax.set_xlim(0,20)
#ax.set_ylim(-7, 0)


majorLocator = MultipleLocator(2)
ax.xaxis.set_major_locator(majorLocator)




#minorLocator = MultipleLocator(0.1)
#minor_loc_2 = MultipleLocator(0.2)
#ax.xaxis.set_minor_locator(minorLocator)
#ax.yaxis.set_minor_locator(minor_loc_2)
#ax.set_xscale('log')
#displaying the figure
plt.show()

fig.savefig('igf_10km.pdf')

