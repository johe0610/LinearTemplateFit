
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                   0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                                                                     [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                                                                     [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                                                                     [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                                                                      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                      1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
                                                                                     [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                      1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                      1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
                                                                                     [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                      1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                      1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
                                                                                     [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                      1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
