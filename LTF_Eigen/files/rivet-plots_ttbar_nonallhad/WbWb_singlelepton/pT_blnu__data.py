
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 1.900000e+02, 2.600000e+02,
                                                                                   3.500000e+02, 4.700000e+02, 6.300000e+02, 8.350000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 1.900000e+02, 2.600000e+02,
                                                                                   3.500000e+02, 4.700000e+02, 6.300000e+02, 8.350000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 1.900000e+02, 2.600000e+02,
                                                                                   3.500000e+02, 4.700000e+02, 6.300000e+02, 8.350000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.200000e+02,
                                                                                   3.000000e+02, 4.000000e+02, 5.400000e+02, 7.200000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.200000e+02,
                                                                                   3.000000e+02, 4.000000e+02, 5.400000e+02, 7.200000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.200000e+02,
                                                                                   3.000000e+02, 4.000000e+02, 5.400000e+02, 7.200000e+02, 9.500000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.097095e-02, 5.322039e-02, 8.770497e-02, 6.297526e-02, 3.124510e-02,
                                                                                   1.154581e-02, 3.109881e-03, 6.377957e-04, 1.110992e-04],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.072874e-02, 5.328498e-02, 8.733870e-02, 6.223414e-02, 3.089857e-02,
                                                                                   1.134242e-02, 3.037197e-03, 6.066214e-04, 1.031586e-04],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.048764e-02, 5.075672e-02, 8.371660e-02, 6.048623e-02, 3.006489e-02,
                                                                                   1.113259e-02, 3.004722e-03, 6.048233e-04, 1.016045e-04],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.113947e-05, 2.454284e-05, 2.879169e-05, 2.444436e-05, 1.495812e-05,
                                                                                      8.173031e-06, 3.617642e-06, 1.464977e-06, 5.460186e-07],
                                                                                     [1.113947e-05, 2.454284e-05, 2.879169e-05, 2.444436e-05, 1.495812e-05,
                                                                                      8.173031e-06, 3.617642e-06, 1.464977e-06, 5.460186e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.151358e-05, 4.795971e-05, 5.610986e-05, 4.746159e-05, 2.905296e-05,
                                                                                      1.582784e-05, 6.984877e-06, 2.797801e-06, 1.030496e-06],
                                                                                     [2.151358e-05, 4.795971e-05, 5.610986e-05, 4.746159e-05, 2.905296e-05,
                                                                                      1.582784e-05, 6.984877e-06, 2.797801e-06, 1.030496e-06],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.041361e-05, 4.492499e-05, 5.272475e-05, 4.489856e-05, 2.750124e-05,
                                                                                      1.504298e-05, 6.663053e-06, 2.682118e-06, 9.844463e-07],
                                                                                     [2.041361e-05, 4.492499e-05, 5.272475e-05, 4.489856e-05, 2.750124e-05,
                                                                                      1.504298e-05, 6.663053e-06, 2.682118e-06, 9.844463e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.779226e-01, 1.001214e+00, 9.958238e-01, 9.882316e-01, 9.889093e-01,
                                                                                   9.823841e-01, 9.766280e-01, 9.511218e-01, 9.285269e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.559464e-01, 9.537082e-01, 9.545252e-01, 9.604761e-01, 9.622274e-01,
                                                                                   9.642104e-01, 9.661855e-01, 9.483026e-01, 9.145385e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.015361e-03, 4.611548e-04, 3.282789e-04, 3.881581e-04, 4.787349e-04,
                                                                                      7.078785e-04, 1.163273e-03, 2.296938e-03, 4.914694e-03],
                                                                                     [1.015361e-03, 4.611548e-04, 3.282789e-04, 3.881581e-04, 4.787349e-04,
                                                                                      7.078785e-04, 1.163273e-03, 2.296938e-03, 4.914694e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.960959e-03, 9.011529e-04, 6.397569e-04, 7.536545e-04, 9.298405e-04,
                                                                                      1.370873e-03, 2.246027e-03, 4.386673e-03, 9.275458e-03],
                                                                                     [1.960959e-03, 9.011529e-04, 6.397569e-04, 7.536545e-04, 9.298405e-04,
                                                                                      1.370873e-03, 2.246027e-03, 4.386673e-03, 9.275458e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.860697e-03, 8.441312e-04, 6.011603e-04, 7.129555e-04, 8.801777e-04,
                                                                                      1.302895e-03, 2.142543e-03, 4.205293e-03, 8.860967e-03],
                                                                                     [1.860697e-03, 8.441312e-04, 6.011603e-04, 7.129555e-04, 8.801777e-04,
                                                                                      1.302895e-03, 2.142543e-03, 4.205293e-03, 8.860967e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
