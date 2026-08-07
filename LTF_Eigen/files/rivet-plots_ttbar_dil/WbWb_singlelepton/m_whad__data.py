
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.800000e+01, 7.600000e+01, 8.400000e+01, 9.200000e+01],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.800000e+01, 7.600000e+01, 8.400000e+01, 9.200000e+01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.800000e+01, 7.600000e+01, 8.400000e+01, 9.200000e+01],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.400000e+01, 7.200000e+01, 8.000000e+01, 8.800000e+01, 9.600000e+01],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.400000e+01, 7.200000e+01, 8.000000e+01, 8.800000e+01, 9.600000e+01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.400000e+01, 7.200000e+01, 8.000000e+01, 8.800000e+01, 9.600000e+01],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.077186e-02, 1.005314e-02, 9.549469e-03, 9.143707e-03],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.097705e-02, 1.032809e-02, 9.830995e-03, 9.360848e-03],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.055670e-02, 9.925024e-03, 9.395727e-03, 9.061180e-03],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.378954e-05, 1.332114e-05, 1.298434e-05, 1.270416e-05],
                                                                                     [1.378954e-05, 1.332114e-05, 1.298434e-05, 1.270416e-05],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.785557e-05, 3.671934e-05, 3.582444e-05, 3.496546e-05],
                                                                                     [3.785557e-05, 3.671934e-05, 3.582444e-05, 3.496546e-05],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.562122e-05, 3.454046e-05, 3.360897e-05, 3.301953e-05],
                                                                                     [3.562122e-05, 3.454046e-05, 3.360897e-05, 3.301953e-05],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.019049e+00, 1.027350e+00, 1.029481e+00, 1.023748e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.800257e-01, 9.872561e-01, 9.839005e-01, 9.909744e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.280145e-03, 1.325073e-03, 1.359692e-03, 1.389388e-03],
                                                                                     [1.280145e-03, 1.325073e-03, 1.359692e-03, 1.389388e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.514302e-03, 3.652524e-03, 3.751459e-03, 3.823992e-03],
                                                                                     [3.514302e-03, 3.652524e-03, 3.751459e-03, 3.823992e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.306877e-03, 3.435788e-03, 3.519460e-03, 3.611175e-03],
                                                                                     [3.306877e-03, 3.435788e-03, 3.519460e-03, 3.611175e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
