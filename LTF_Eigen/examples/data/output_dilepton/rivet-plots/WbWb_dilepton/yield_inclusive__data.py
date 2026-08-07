
import numpy as np
from numpy import nan

add_legend_handle = [
    'S3beta_Cluster_WbWb_dilepton_170.yoda.gz',
    'S3beta_Cluster_WbWb_dilepton_160.yoda.gz',
    'S3beta_Cluster_WbWb_dilepton_165.yoda.gz',
    'S3beta_Cluster_WbWb_dilepton_175.yoda.gz',
    'S3beta_Cluster_WbWb_dilepton_180.yoda.gz'
]

xpoints = {
    'S3beta_Cluster_WbWb_dilepton_170.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_160.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_165.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_175.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_180.yoda.gz' : [1.0],
}
xedges = {
    'S3beta_Cluster_WbWb_dilepton_170.yoda.gz' : [0.5, 1.5],
    'S3beta_Cluster_WbWb_dilepton_160.yoda.gz' : [0.5, 1.5],
    'S3beta_Cluster_WbWb_dilepton_165.yoda.gz' : [0.5, 1.5],
    'S3beta_Cluster_WbWb_dilepton_175.yoda.gz' : [0.5, 1.5],
    'S3beta_Cluster_WbWb_dilepton_180.yoda.gz' : [0.5, 1.5],
}
ref_xerrs = [
  [abs(xpoints['S3beta_Cluster_WbWb_dilepton_170.yoda.gz'][i]   - xedges['S3beta_Cluster_WbWb_dilepton_170.yoda.gz'][i]) for i in range(len(xpoints['S3beta_Cluster_WbWb_dilepton_170.yoda.gz']))],
  [abs(xedges['S3beta_Cluster_WbWb_dilepton_170.yoda.gz'][i+1] - xpoints['S3beta_Cluster_WbWb_dilepton_170.yoda.gz'][i]) for i in range(len(xpoints['S3beta_Cluster_WbWb_dilepton_170.yoda.gz']))]
]

yvals = {
    'S3beta_Cluster_WbWb_dilepton_170.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_160.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_165.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_175.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_180.yoda.gz' : [1.0],
}
xerrs = {
    'S3beta_Cluster_WbWb_dilepton_170.yoda.gz' : [
        [0.5],
        [0.5],
    ],
    'S3beta_Cluster_WbWb_dilepton_160.yoda.gz' : [
        [0.5],
        [0.5],
    ],
    'S3beta_Cluster_WbWb_dilepton_165.yoda.gz' : [
        [0.5],
        [0.5],
    ],
    'S3beta_Cluster_WbWb_dilepton_175.yoda.gz' : [
        [0.5],
        [0.5],
    ],
    'S3beta_Cluster_WbWb_dilepton_180.yoda.gz' : [
        [0.5],
        [0.5],
    ],
}
yerrs = {
    'S3beta_Cluster_WbWb_dilepton_170.yoda.gz' : [
        [0.007340600111707489],
        [0.007340600111707489],
    ],
    'S3beta_Cluster_WbWb_dilepton_160.yoda.gz' : [
        [0.007534637350264444],
        [0.007534637350264444],
    ],
    'S3beta_Cluster_WbWb_dilepton_165.yoda.gz' : [
        [0.006694051090333865],
        [0.006694051090333865],
    ],
    'S3beta_Cluster_WbWb_dilepton_175.yoda.gz' : [
        [0.00778006940843075],
        [0.00778006940843075],
    ],
    'S3beta_Cluster_WbWb_dilepton_180.yoda.gz' : [
        [0.007049213431298559],
        [0.007049213431298559],
    ],
}
variation_yvals = {
}


# lists for ratio plot
ratio0_yvals = {
    'S3beta_Cluster_WbWb_dilepton_170.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_160.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_165.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_175.yoda.gz' : [1.0],
    'S3beta_Cluster_WbWb_dilepton_180.yoda.gz' : [1.0],
}
ratio0_yerrs = {
    'S3beta_Cluster_WbWb_dilepton_170.yoda.gz' : [
        [0.007340600111707489],
        [0.007340600111707489],
    ],
    'S3beta_Cluster_WbWb_dilepton_160.yoda.gz' : [
        [0.007534637350264444],
        [0.007534637350264444],
    ],
    'S3beta_Cluster_WbWb_dilepton_165.yoda.gz' : [
        [0.006694051090333865],
        [0.006694051090333865],
    ],
    'S3beta_Cluster_WbWb_dilepton_175.yoda.gz' : [
        [0.00778006940843075],
        [0.00778006940843075],
    ],
    'S3beta_Cluster_WbWb_dilepton_180.yoda.gz' : [
        [0.007049213431298559],
        [0.007049213431298559],
    ],
}
ratio0_variation_vals = {
}
ratio_band_edges = {
}