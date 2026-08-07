
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 2.000000e+02, 2.800000e+02,
                                                                                   3.800000e+02, 5.050000e+02, 7.350000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 2.000000e+02, 2.800000e+02,
                                                                                   3.800000e+02, 5.050000e+02, 7.350000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 2.000000e+02, 2.800000e+02,
                                                                                   3.800000e+02, 5.050000e+02, 7.350000e+02, 9.500000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.400000e+02,
                                                                                   3.200000e+02, 4.400000e+02, 5.700000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.400000e+02,
                                                                                   3.200000e+02, 4.400000e+02, 5.700000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.400000e+02,
                                                                                   3.200000e+02, 4.400000e+02, 5.700000e+02, 9.000000e+02, 1.000000e+03],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.535453e-04, 1.050220e-03, 2.778553e-03, 1.622904e-03, 3.611462e-04,
                                                                                   7.730686e-05, 1.492205e-05, 1.883645e-06, 2.268580e-07],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.551366e-04, 1.099011e-03, 2.879979e-03, 1.635146e-03, 3.616890e-04,
                                                                                   7.839739e-05, 1.496198e-05, 1.880084e-06, 1.760037e-07],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.484403e-04, 1.020626e-03, 2.700064e-03, 1.627894e-03, 3.617751e-04,
                                                                                   7.793822e-05, 1.511366e-05, 1.899689e-06, 2.099112e-07],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.581515e-07, 1.720650e-06, 2.554901e-06, 1.692343e-06, 8.009330e-07,
                                                                                      3.040836e-07, 1.297402e-07, 2.942047e-08, 1.868354e-08],
                                                                                     [6.581515e-07, 1.720650e-06, 2.554901e-06, 1.692343e-06, 8.009330e-07,
                                                                                      3.040836e-07, 1.297402e-07, 2.942047e-08, 1.868354e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.801280e-06, 4.786504e-06, 7.073770e-06, 4.618682e-06, 2.181798e-06,
                                                                                      8.337555e-07, 3.539715e-07, 8.130406e-08, 5.380449e-08],
                                                                                     [1.801280e-06, 4.786504e-06, 7.073770e-06, 4.618682e-06, 2.181798e-06,
                                                                                      8.337555e-07, 3.539715e-07, 8.130406e-08, 5.380449e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.687883e-06, 4.425701e-06, 6.573151e-06, 4.423368e-06, 2.091850e-06,
                                                                                      7.966748e-07, 3.410683e-07, 7.753374e-08, 4.865669e-08],
                                                                                     [1.687883e-06, 4.425701e-06, 6.573151e-06, 4.423368e-06, 2.091850e-06,
                                                                                      7.966748e-07, 3.410683e-07, 7.753374e-08, 4.865669e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.010364e+00, 1.046458e+00, 1.036503e+00, 1.007543e+00, 1.001503e+00,
                                                                                   1.014107e+00, 1.002676e+00, 9.981095e-01, 7.758320e-01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.667525e-01, 9.718211e-01, 9.717518e-01, 1.003075e+00, 1.001741e+00,
                                                                                   1.008167e+00, 1.012841e+00, 1.008518e+00, 9.252978e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.286367e-03, 1.638371e-03, 9.195077e-04, 1.042787e-03, 2.217753e-03,
                                                                                      3.933462e-03, 8.694529e-03, 1.561890e-02, 8.235786e-02],
                                                                                     [4.286367e-03, 1.638371e-03, 9.195077e-04, 1.042787e-03, 2.217753e-03,
                                                                                      3.933462e-03, 8.694529e-03, 1.561890e-02, 8.235786e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.173126e-02, 4.557620e-03, 2.545847e-03, 2.845937e-03, 6.041315e-03,
                                                                                      1.078501e-02, 2.372137e-02, 4.316315e-02, 2.371725e-01],
                                                                                     [1.173126e-02, 4.557620e-03, 2.545847e-03, 2.845937e-03, 6.041315e-03,
                                                                                      1.078501e-02, 2.372137e-02, 4.316315e-02, 2.371725e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.099274e-02, 4.214070e-03, 2.365674e-03, 2.725588e-03, 5.792253e-03,
                                                                                      1.030536e-02, 2.285667e-02, 4.116155e-02, 2.144808e-01],
                                                                                     [1.099274e-02, 4.214070e-03, 2.365674e-03, 2.725588e-03, 5.792253e-03,
                                                                                      1.030536e-02, 2.285667e-02, 4.116155e-02, 2.144808e-01],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
