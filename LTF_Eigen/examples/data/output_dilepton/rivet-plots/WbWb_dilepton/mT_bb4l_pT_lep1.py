#! /usr/bin/env python

# This Python script was auto-generated using YODA v2.0.1.
# Analysis object: /WbWb_dilepton/mT_bb4l_pT_lep1
# Timestamp: 27-02-2025 (11:02:16)

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg') # comment out for interactive use
import os
import numpy as np

plotDir = os.path.split(os.path.realpath(__file__))[0]
if 'YODA_USER_PLOT_PATH' in globals():
    plot_outdir = globals()['YODA_USER_PLOT_PATH']
else:
    plot_outdir = plotDir

#plot style
plt.style.use(os.path.join(plotDir, '../default.mplstyle'))

plt.rcParams['text.usetex'] = False
plt.style.use(['classic'])
plt.rc('font', family='DejaVu Sans')
# plot metadata
figW, figH = plt.rcParams['figure.figsize']
ax_xLabel = r''
ax_yLabel = r''
ax_zLabel = r''
ax_title  = r''
ax_xScale = 'linear'
ax_yScale = 'log'
ax_zScale = 'linear'
yLims = (10, 100)

# TeX-friendly labels for the legend
labels = [ r"S3beta\_Cluster\_WbWb\_dilepton\_170", r"S3beta\_Cluster\_WbWb\_dilepton\_160", r"S3beta\_Cluster\_WbWb\_dilepton\_165", r"S3beta\_Cluster\_WbWb\_dilepton\_175", r"S3beta\_Cluster\_WbWb\_dilepton\_180" ]

# Adjust canvas width and height
canvasW = 10
canvasH = 10
figW *= canvasW/10.
figH *= canvasH/9.

# Create figure and axis objects
fig, ax = plt.subplots(1, 1)

# Set figure margins
plt.subplots_adjust(
    left   = 1.0 * plt.rcParams['figure.subplot.left'],
    right  = 1.0 * plt.rcParams['figure.subplot.right'],
    top    = 1.0 * plt.rcParams['figure.subplot.top'],
    bottom = 1.0 * plt.rcParams['figure.subplot.bottom'])

ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(
base=10.0, subs=[i for i in np.arange(0, 1, 0.1)], numticks=np.inf))

# Color map: curve index -> color
colors = {0: '#EE3311', 1: '#3366FF', 2: '#109618', 3: '#FF9900', 4: '#990099'}


# the numerical data is stored in a separate file
dataf = dict()
exec(open(os.path.split(__file__)[0] + '/mT_bb4l_pT_lep1__data.py').read(), dataf)

legend_handles = dict() # keep track of handles for the legend
# reference data in main panel
cmap = 'cividis'
xbin_cent = np.unique(dataf['xpoints'])
ybin_cent = np.unique(dataf['ypoints'])
X, Y = np.meshgrid(xbin_cent, ybin_cent)
Z = np.array(dataf['zpoints']).reshape(X.shape[::-1])
pc = ax.pcolormesh(X, Y, Z.T, cmap=cmap, shading='auto')
cbar = fig.colorbar(pc, orientation='vertical', ax=ax, pad=0.01)
cbar.set_label(ax_zLabel)
# style options for curves
# starts at zorder>=5 to draw curve on top of legend

legend_pos = (0.5, 0.97)
ax.legend(legend_handles, labels, loc='upper left', bbox_to_anchor=legend_pos)

ax.set_xlabel(ax_xLabel)
ax.set_ylabel(ax_yLabel, loc='top')
ax.set_title(ax_title, loc='left')
ax.set_xscale(ax_xScale)
ax.set_yscale(ax_yScale)

# tick formatting
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.savefig(os.path.join(plot_outdir, 'mT_bb4l_pT_lep1.pdf'), format='PDF')
plt.savefig(os.path.join(plot_outdir, 'mT_bb4l_pT_lep1.png'), format='PNG')

plt.close(fig)