
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 1.200000e+02, 2.000000e+02, 2.800000e+02, 3.800000e+02,
                                                                                   5.200000e+02, 7.500000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 1.200000e+02, 2.000000e+02, 2.800000e+02, 3.800000e+02,
                                                                                   5.200000e+02, 7.500000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 1.200000e+02, 2.000000e+02, 2.800000e+02, 3.800000e+02,
                                                                                   5.200000e+02, 7.500000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 8.000000e+01, 1.600000e+02, 2.400000e+02, 3.200000e+02,
                                                                                   4.400000e+02, 6.000000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 8.000000e+01, 1.600000e+02, 2.400000e+02, 3.200000e+02,
                                                                                   4.400000e+02, 6.000000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 8.000000e+01, 1.600000e+02, 2.400000e+02, 3.200000e+02,
                                                                                   4.400000e+02, 6.000000e+02, 9.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.791128e-04, 2.046359e-03, 1.365133e-03, 5.222660e-04, 1.823279e-04,
                                                                                   6.035107e-05, 1.396631e-05],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.885980e-04, 2.104880e-03, 1.398074e-03, 5.323821e-04, 1.880073e-04,
                                                                                   6.043103e-05, 1.450772e-05],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.703096e-04, 2.004311e-03, 1.348333e-03, 5.162408e-04, 1.826164e-04,
                                                                                   6.056383e-05, 1.371261e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.009728e-06, 1.898897e-06, 1.553229e-06, 9.619577e-07, 4.643981e-07,
                                                                                      2.311934e-07, 8.106140e-08],
                                                                                     [1.009728e-06, 1.898897e-06, 1.553229e-06, 9.619577e-07, 4.643981e-07,
                                                                                      2.311934e-07, 8.106140e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.768730e-06, 5.237337e-06, 4.274803e-06, 2.641450e-06, 1.283197e-06,
                                                                                      6.296058e-07, 2.247087e-07],
                                                                                     [2.768730e-06, 5.237337e-06, 4.274803e-06, 2.641450e-06, 1.283197e-06,
                                                                                      6.296058e-07, 2.247087e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.614774e-06, 4.905026e-06, 4.028394e-06, 2.495647e-06, 1.211649e-06,
                                                                                      6.053035e-07, 2.098148e-07],
                                                                                     [2.614774e-06, 4.905026e-06, 4.028394e-06, 2.495647e-06, 1.211649e-06,
                                                                                      6.053035e-07, 2.098148e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.016379e+00, 1.028598e+00, 1.024130e+00, 1.019370e+00, 1.031149e+00,
                                                                                   1.001325e+00, 1.038765e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.847988e-01, 9.794523e-01, 9.876935e-01, 9.884634e-01, 1.001582e+00,
                                                                                   1.003525e+00, 9.818349e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.743577e-03, 9.279393e-04, 1.137786e-03, 1.841892e-03, 2.547049e-03,
                                                                                      3.830809e-03, 5.804067e-03],
                                                                                     [1.743577e-03, 9.279393e-04, 1.137786e-03, 1.841892e-03, 2.547049e-03,
                                                                                      3.830809e-03, 5.804067e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.780986e-03, 2.559344e-03, 3.131419e-03, 5.057672e-03, 7.037853e-03,
                                                                                      1.043239e-02, 1.608934e-02],
                                                                                     [4.780986e-03, 2.559344e-03, 3.131419e-03, 5.057672e-03, 7.037853e-03,
                                                                                      1.043239e-02, 1.608934e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.515138e-03, 2.396953e-03, 2.950917e-03, 4.778498e-03, 6.645439e-03,
                                                                                      1.002971e-02, 1.502292e-02],
                                                                                     [4.515138e-03, 2.396953e-03, 2.950917e-03, 4.778498e-03, 6.645439e-03,
                                                                                      1.002971e-02, 1.502292e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
