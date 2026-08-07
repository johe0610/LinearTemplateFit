
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.000000e+01, 9.000000e+01, 1.300000e+02, 1.700000e+02, 2.100000e+02,
                                                                                   2.600000e+02, 3.200000e+02, 4.000000e+02, 5.250000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.000000e+01, 9.000000e+01, 1.300000e+02, 1.700000e+02, 2.100000e+02,
                                                                                   2.600000e+02, 3.200000e+02, 4.000000e+02, 5.250000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.000000e+01, 9.000000e+01, 1.300000e+02, 1.700000e+02, 2.100000e+02,
                                                                                   2.600000e+02, 3.200000e+02, 4.000000e+02, 5.250000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 7.000000e+01, 1.100000e+02, 1.500000e+02, 1.900000e+02,
                                                                                   2.300000e+02, 2.900000e+02, 3.500000e+02, 4.500000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 7.000000e+01, 1.100000e+02, 1.500000e+02, 1.900000e+02,
                                                                                   2.300000e+02, 2.900000e+02, 3.500000e+02, 4.500000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 7.000000e+01, 1.100000e+02, 1.500000e+02, 1.900000e+02,
                                                                                   2.300000e+02, 2.900000e+02, 3.500000e+02, 4.500000e+02, 6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.981937e-03, 2.397993e-03, 9.135441e-04, 3.609470e-04, 1.491944e-04,
                                                                                   5.638411e-05, 1.951355e-05, 5.837864e-06, 1.227610e-06],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.144030e-03, 2.457347e-03, 9.321355e-04, 3.618129e-04, 1.507803e-04,
                                                                                   5.741605e-05, 1.931692e-05, 5.903446e-06, 1.298216e-06],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.869952e-03, 2.375952e-03, 9.027999e-04, 3.592167e-04, 1.502700e-04,
                                                                                   5.534595e-05, 1.946647e-05, 5.782602e-06, 1.208564e-06],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.591295e-06, 2.910277e-06, 1.798938e-06, 1.132369e-06, 7.300950e-07,
                                                                                      3.674727e-07, 2.166982e-07, 9.268474e-08, 3.469244e-08],
                                                                                     [4.591295e-06, 2.910277e-06, 1.798938e-06, 1.132369e-06, 7.300950e-07,
                                                                                      3.674727e-07, 2.166982e-07, 9.268474e-08, 3.469244e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.265413e-05, 8.012280e-06, 4.941911e-06, 3.088744e-06, 1.993838e-06,
                                                                                      1.005396e-06, 5.927373e-07, 2.553775e-07, 9.716969e-08],
                                                                                     [1.265413e-05, 8.012280e-06, 4.941911e-06, 3.088744e-06, 1.993838e-06,
                                                                                      1.005396e-06, 5.927373e-07, 2.553775e-07, 9.716969e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.186937e-05, 7.562058e-06, 4.665143e-06, 2.949225e-06, 1.908346e-06,
                                                                                      9.491921e-07, 5.672180e-07, 2.391292e-07, 8.950404e-08],
                                                                                     [1.186937e-05, 7.562058e-06, 4.665143e-06, 2.949225e-06, 1.908346e-06,
                                                                                      9.491921e-07, 5.672180e-07, 2.391292e-07, 8.950404e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.027097e+00, 1.024752e+00, 1.020351e+00, 1.002399e+00, 1.010630e+00,
                                                                                   1.018302e+00, 9.899234e-01, 1.011234e+00, 1.057515e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.812795e-01, 9.908086e-01, 9.882390e-01, 9.952062e-01, 1.007209e+00,
                                                                                   9.815877e-01, 9.975873e-01, 9.905339e-01, 9.844853e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.675265e-04, 1.213630e-03, 1.969186e-03, 3.137217e-03, 4.893582e-03,
                                                                                      6.517310e-03, 1.110501e-02, 1.587648e-02, 2.826015e-02],
                                                                                     [7.675265e-04, 1.213630e-03, 1.969186e-03, 3.137217e-03, 4.893582e-03,
                                                                                      6.517310e-03, 1.110501e-02, 1.587648e-02, 2.826015e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.115390e-03, 3.341244e-03, 5.409603e-03, 8.557334e-03, 1.336403e-02,
                                                                                      1.783119e-02, 3.037568e-02, 4.374502e-02, 7.915355e-02],
                                                                                     [2.115390e-03, 3.341244e-03, 5.409603e-03, 8.557334e-03, 1.336403e-02,
                                                                                      1.783119e-02, 3.037568e-02, 4.374502e-02, 7.915355e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.984202e-03, 3.153495e-03, 5.106642e-03, 8.170798e-03, 1.279100e-02,
                                                                                      1.683439e-02, 2.906790e-02, 4.096176e-02, 7.290918e-02],
                                                                                     [1.984202e-03, 3.153495e-03, 5.106642e-03, 8.170798e-03, 1.279100e-02,
                                                                                      1.683439e-02, 2.906790e-02, 4.096176e-02, 7.290918e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
