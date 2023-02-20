#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Jasmine Hansen, 2023
"""


import pandas as pd
import os 
import matplotlib.pyplot as plt
import numpy as np
import itertools

# read in data 
comb_data = pd.read_csv('comb_results.csv')

tomo = comb_data[comb_data['site_id'] == 'TOMO']
berp = comb_data[comb_data['site_id'] == 'BERP']
mrtp = comb_data[comb_data['site_id'] == 'MRTP']
mtak = comb_data[comb_data['site_id'] == 'MTAK']
inmn = comb_data[comb_data['site_id'] == 'INMN']
thur = comb_data[comb_data['site_id'] == 'THUR']
lply = comb_data[comb_data['site_id'] == 'LPLY']
sltr = comb_data[comb_data['site_id'] == 'SLTR']

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

names = ['0.32', '0.86', '2.27', '6.00']
#fig, ax0 = plt.subplots(1,1)


#ax0.errorbar(tomo['run_val'][3], tomo['V55_mean'][3], yerr = tomo['v55_std'][3], fmt='o', marker='s', c=CB_color_cycle[8], markeredgecolor = CB_color_cycle[0], ecolor = CB_color_cycle[0], capsize=(2), lw=1, markersize=7)



fig, ((ax0, ax2, ax3, ax5), (ax1, ax4, ax7, ax6)) = plt.subplots(nrows=2, ncols=4, 
                                    figsize=(12, 6))

#fig, ((ax0, ax2), (ax3, ax5), (ax1, ax4), (ax7, ax6)) = plt.subplots(nrows=4, ncols=2, 
#                                    figsize=(8.5, 11))

fig.tight_layout(h_pad=2)

ax0.errorbar(tomo['run_val'][3], tomo['V55_mean'][3], yerr = tomo['v55_std'][3], fmt='o', marker='s', c=CB_color_cycle[0], markeredgecolor = CB_color_cycle[0], ecolor = CB_color_cycle[0], capsize=(2), lw=1, markersize=7)
ax0.errorbar(tomo['run_val'][2], tomo['V55_mean'][2], yerr = tomo['v55_std'][2], fmt='o', marker='D', c=CB_color_cycle[0], markeredgecolor = CB_color_cycle[0], ecolor = CB_color_cycle[0], capsize=(2), lw=1, markersize=7)
ax0.errorbar(tomo['run_val'][1], tomo['V55_mean'][1], yerr = tomo['v55_std'][1], fmt='o', marker='^', c=CB_color_cycle[0], markeredgecolor = CB_color_cycle[0], ecolor = CB_color_cycle[0], capsize=(2), lw=1, markersize=7)
ax0.errorbar(tomo['run_val'][0], tomo['V55_mean'][0], yerr = tomo['v55_std'][0], fmt='o', marker='o', c=CB_color_cycle[0], markeredgecolor = CB_color_cycle[0], ecolor = CB_color_cycle[0], capsize=(2), lw=1, markersize=7)

ax1.errorbar(berp['run_val'][7], berp['V55_mean'][7], yerr = berp['v55_std'][7], fmt='o', marker='s', c=CB_color_cycle[1], markeredgecolor = CB_color_cycle[1], ecolor = CB_color_cycle[1], capsize=(2), lw=1, markersize=7)
ax1.errorbar(berp['run_val'][6], berp['V55_mean'][6], yerr = berp['v55_std'][6], fmt='o', marker='D', c=CB_color_cycle[1], markeredgecolor = CB_color_cycle[1], ecolor = CB_color_cycle[1], capsize=(2), lw=1, markersize=7)
ax1.errorbar(berp['run_val'][5], berp['V55_mean'][5], yerr = berp['v55_std'][5], fmt='o', marker='^', c=CB_color_cycle[1], markeredgecolor = CB_color_cycle[1], ecolor = CB_color_cycle[1], capsize=(2), lw=1, markersize=7)
ax1.errorbar(berp['run_val'][4], berp['V55_mean'][4], yerr = berp['v55_std'][4], fmt='o', marker='o', c=CB_color_cycle[1], markeredgecolor = CB_color_cycle[1], ecolor = CB_color_cycle[1], capsize=(2), lw=1, markersize=7)

ax2.errorbar(sltr['run_val'][11], sltr['V55_mean'][11], yerr = sltr['v55_std'][11], fmt='o', marker='s', c=CB_color_cycle[2], markeredgecolor = CB_color_cycle[2], ecolor = CB_color_cycle[2], capsize=(2), lw=1, markersize=7)
ax2.errorbar(sltr['run_val'][10], sltr['V55_mean'][10], yerr = sltr['v55_std'][10], fmt='o', marker='D', c=CB_color_cycle[2], markeredgecolor = CB_color_cycle[2], ecolor = CB_color_cycle[2], capsize=(2), lw=1, markersize=7)
ax2.errorbar(sltr['run_val'][9], sltr['V55_mean'][9], yerr = sltr['v55_std'][9], fmt='o', marker='^', c=CB_color_cycle[2], markeredgecolor = CB_color_cycle[2], ecolor = CB_color_cycle[2], capsize=(2), lw=1, markersize=7)
ax2.errorbar(sltr['run_val'][8], sltr['V55_mean'][8], yerr = sltr['v55_std'][8], fmt='o', marker='o', c=CB_color_cycle[2], markeredgecolor = CB_color_cycle[2], ecolor = CB_color_cycle[2], capsize=(2), lw=1, markersize=7)

ax3.errorbar(mtak['run_val'][15], mtak['V55_mean'][15], yerr = mtak['v55_std'][15], fmt='o', marker='s', c=CB_color_cycle[5], markeredgecolor = CB_color_cycle[5], ecolor = CB_color_cycle[5], capsize=(2), lw=1, markersize=7)
ax3.errorbar(mtak['run_val'][14], mtak['V55_mean'][14], yerr = mtak['v55_std'][14], fmt='o', marker='D', c=CB_color_cycle[5], markeredgecolor = CB_color_cycle[5], ecolor = CB_color_cycle[5], capsize=(2), lw=1, markersize=7)
ax3.errorbar(mtak['run_val'][13], mtak['V55_mean'][13], yerr = mtak['v55_std'][13], fmt='o', marker='^', c=CB_color_cycle[5], markeredgecolor = CB_color_cycle[5], ecolor = CB_color_cycle[5], capsize=(2), lw=1, markersize=7)
ax3.errorbar(mtak['run_val'][12], mtak['V55_mean'][12], yerr = mtak['v55_std'][12], fmt='o', marker='o', c=CB_color_cycle[5], markeredgecolor = CB_color_cycle[5], ecolor = CB_color_cycle[5], capsize=(2), lw=1, markersize=7)

ax4.errorbar(mrtp['run_val'][19], mrtp['V55_mean'][19], yerr = mrtp['v55_std'][19], fmt='o', marker='s', c='brown', markeredgecolor = 'brown', ecolor = 'brown', capsize=(2), lw=1, markersize=7)
ax4.errorbar(mrtp['run_val'][18], mrtp['V55_mean'][18], yerr = mrtp['v55_std'][18], fmt='o', marker='D', c='brown', markeredgecolor = 'brown', ecolor = 'brown', capsize=(2), lw=1, markersize=7)
ax4.errorbar(mrtp['run_val'][17], mrtp['V55_mean'][17], yerr = mrtp['v55_std'][17], fmt='o', marker='^', c='brown', markeredgecolor = 'brown', ecolor = 'brown', capsize=(2), lw=1, markersize=7)
ax4.errorbar(mrtp['run_val'][16], mrtp['V55_mean'][16], yerr = mrtp['v55_std'][16], fmt='o', marker='o', c='brown', markeredgecolor = 'brown', ecolor = 'brown', capsize=(2), lw=1, markersize=7)

ax5.errorbar(inmn['run_val'][23], inmn['V55_mean'][23], yerr = inmn['v55_std'][23], fmt='o', marker='s', c='teal', markeredgecolor = 'teal', ecolor = 'teal', capsize=(2), lw=1, markersize=7)
ax5.errorbar(inmn['run_val'][22], inmn['V55_mean'][22], yerr = inmn['v55_std'][22], fmt='o', marker='D', c='teal', markeredgecolor = 'teal', ecolor = 'teal', capsize=(2), lw=1, markersize=7)
ax5.errorbar(inmn['run_val'][21], inmn['V55_mean'][21], yerr = inmn['v55_std'][21], fmt='o', marker='^', c='teal', markeredgecolor = 'teal', ecolor = 'teal', capsize=(2), lw=1, markersize=7)
ax5.errorbar(inmn['run_val'][20], inmn['V55_mean'][20], yerr = inmn['v55_std'][20], fmt='o', marker='o', c='teal', markeredgecolor = 'teal', ecolor = 'teal', capsize=(2), lw=1, markersize=7)

ax6.errorbar(thur['run_val'][27], thur['V55_mean'][27], yerr = thur['v55_std'][27], fmt='o', marker='s', c='k', markeredgecolor = 'k', ecolor = 'k', capsize=(2), lw=1, markersize=7)
ax6.errorbar(thur['run_val'][26], thur['V55_mean'][26], yerr = thur['v55_std'][26], fmt='o', marker='D', c='k', markeredgecolor = 'k', ecolor = 'k', capsize=(2), lw=1, markersize=7)
ax6.errorbar(thur['run_val'][25], thur['V55_mean'][25], yerr = thur['v55_std'][25], fmt='o', marker='^', c='k', markeredgecolor = 'k', ecolor = 'k', capsize=(2), lw=1, markersize=7)
ax6.errorbar(thur['run_val'][24], thur['V55_mean'][24], yerr = thur['v55_std'][24], fmt='o', marker='o', c='k', markeredgecolor = 'k', ecolor = 'k', capsize=(2), lw=1, markersize=7)

ax7.errorbar(lply['run_val'][31], lply['V55_mean'][31], yerr = lply['v55_std'][31], fmt='o', marker='s', c='r', markeredgecolor = 'r', ecolor = 'r', capsize=(2), lw=1, markersize=7)
ax7.errorbar(lply['run_val'][30], lply['V55_mean'][30], yerr = lply['v55_std'][30], fmt='o', marker='D', c='r', markeredgecolor = 'r', ecolor = 'r', capsize=(2), lw=1, markersize=7)
ax7.errorbar(lply['run_val'][29], lply['V55_mean'][29], yerr = lply['v55_std'][29], fmt='o', marker='^', c='r', markeredgecolor = 'r', ecolor = 'r', capsize=(2), lw=1, markersize=7)
ax7.errorbar(lply['run_val'][28], lply['V55_mean'][28], yerr = lply['v55_std'][28], fmt='o', marker='o', c='r', markeredgecolor = 'r', ecolor = 'r', capsize=(2), lw=1, markersize=7)

ax0.set_title('TOMO')
ax1.set_title('BERP')
ax2.set_title('SLTR')
ax3.set_title('MTAK')
ax4.set_title('MRTP')
ax5.set_title('INMN')
ax6.set_title('THUR')
ax7.set_title('LPLY')

ax0.set_ylim(10.3,11.3)
ax2.set_ylim(14.6,15.7)
ax3.set_ylim(7.5,8.5)
ax5.set_ylim(7.5,8.5)

ax1.set_ylim(6.3,7.3)
ax4.set_ylim(4.1,5.1)
ax7.set_ylim(2.5,3.5)
ax6.set_ylim(1.5,2.5)



ax0_ticks = np.arange(10,12,1)
#ax0.set_yticks(ax0_ticks)

fig.savefig('elastic_subplots_set_scale_3.pdf')


