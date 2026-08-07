
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.700000e+02, 2.700000e+02,
                                                                                   4.700000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.700000e+02, 2.700000e+02,
                                                                                   4.700000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.700000e+02, 2.700000e+02,
                                                                                   4.700000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 2.000000e+02,
                                                                                   3.400000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 2.000000e+02,
                                                                                   3.400000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 2.000000e+02,
                                                                                   3.400000e+02, 6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.526725e-03, 2.876889e-03, 8.712613e-04, 2.362669e-04, 3.249645e-05,
                                                                                   1.564003e-06],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.745577e-03, 2.937014e-03, 8.704262e-04, 2.329447e-04, 3.195972e-05,
                                                                                   1.553638e-06],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.362425e-03, 2.865492e-03, 8.714848e-04, 2.391012e-04, 3.335662e-05,
                                                                                   1.607504e-06],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.126907e-06, 3.186300e-06, 1.757939e-06, 7.499826e-07, 1.833589e-07,
                                                                                      3.028478e-08],
                                                                                     [5.126907e-06, 3.186300e-06, 1.757939e-06, 7.499826e-07, 1.833589e-07,
                                                                                      3.028478e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.417493e-05, 8.754819e-06, 4.779458e-06, 2.027341e-06, 4.959485e-07,
                                                                                      8.260584e-08],
                                                                                     [1.417493e-05, 8.754819e-06, 4.779458e-06, 2.027341e-06, 4.959485e-07,
                                                                                      8.260584e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.321059e-05, 8.300020e-06, 4.585881e-06, 1.967861e-06, 4.844522e-07,
                                                                                      8.140661e-08],
                                                                                     [1.321059e-05, 8.300020e-06, 4.585881e-06, 1.967861e-06, 4.844522e-07,
                                                                                      8.140661e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.033532e+00, 1.020899e+00, 9.990415e-01, 9.859388e-01, 9.834834e-01,
                                                                                   9.933728e-01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.748266e-01, 9.960384e-01, 1.000257e+00, 1.011996e+00, 1.026470e+00,
                                                                                   1.027814e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.855252e-04, 1.107551e-03, 2.017694e-03, 3.174302e-03, 5.642429e-03,
                                                                                      1.936363e-02],
                                                                                     [7.855252e-04, 1.107551e-03, 2.017694e-03, 3.174302e-03, 5.642429e-03,
                                                                                      1.936363e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.171829e-03, 3.043155e-03, 5.485677e-03, 8.580724e-03, 1.526162e-02,
                                                                                      5.281693e-02],
                                                                                     [2.171829e-03, 3.043155e-03, 5.485677e-03, 8.580724e-03, 1.526162e-02,
                                                                                      5.281693e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.024076e-03, 2.885068e-03, 5.263497e-03, 8.328975e-03, 1.490785e-02,
                                                                                      5.205016e-02],
                                                                                     [2.024076e-03, 2.885068e-03, 5.263497e-03, 8.328975e-03, 1.490785e-02,
                                                                                      5.205016e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
