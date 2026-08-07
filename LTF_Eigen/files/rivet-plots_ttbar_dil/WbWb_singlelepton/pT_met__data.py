
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 8.500000e+01, 1.650000e+02, 2.650000e+02, 4.100000e+02,
                                                                                   7.000000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 8.500000e+01, 1.650000e+02, 2.650000e+02, 4.100000e+02,
                                                                                   7.000000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 8.500000e+01, 1.650000e+02, 2.650000e+02, 4.100000e+02,
                                                                                   7.000000e+02, 9.500000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.200000e+02, 2.100000e+02, 3.200000e+02,
                                                                                   5.000000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.200000e+02, 2.100000e+02, 3.200000e+02,
                                                                                   5.000000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.200000e+02, 2.100000e+02, 3.200000e+02,
                                                                                   5.000000e+02, 9.000000e+02, 1.000000e+03],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.680855e-03, 3.001266e-03, 9.035218e-04, 1.523496e-04, 2.453819e-05,
                                                                                   2.136194e-06, 1.680584e-07],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.723262e-03, 3.077938e-03, 9.236707e-04, 1.560964e-04, 2.458822e-05,
                                                                                   2.237240e-06, 1.967925e-07],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.638669e-03, 2.954456e-03, 9.009163e-04, 1.496860e-04, 2.453976e-05,
                                                                                   2.175636e-06, 1.717975e-07],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.177680e-06, 2.459905e-06, 1.191626e-06, 4.431602e-07, 1.391079e-07,
                                                                                      2.746292e-08, 1.534157e-08],
                                                                                     [2.177680e-06, 2.459905e-06, 1.191626e-06, 4.431602e-07, 1.391079e-07,
                                                                                      2.746292e-08, 1.534157e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.995704e-06, 6.776218e-06, 3.275780e-06, 1.220220e-06, 3.792172e-07,
                                                                                      7.646406e-08, 4.514730e-08],
                                                                                     [5.995704e-06, 6.776218e-06, 3.275780e-06, 1.220220e-06, 3.792172e-07,
                                                                                      7.646406e-08, 4.514730e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.613075e-06, 6.369078e-06, 3.104145e-06, 1.147620e-06, 3.630926e-07,
                                                                                      7.243675e-08, 4.049308e-08],
                                                                                     [5.613075e-06, 6.369078e-06, 3.104145e-06, 1.147620e-06, 3.630926e-07,
                                                                                      7.243675e-08, 4.049308e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.025229e+00, 1.025547e+00, 1.022300e+00, 1.024593e+00, 1.002039e+00,
                                                                                   1.047302e+00, 1.170977e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.749021e-01, 9.844032e-01, 9.971163e-01, 9.825165e-01, 1.000064e+00,
                                                                                   1.018464e+00, 1.022249e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.295579e-03, 8.196225e-04, 1.318868e-03, 2.908837e-03, 5.669037e-03,
                                                                                      1.285600e-02, 9.128714e-02],
                                                                                     [1.295579e-03, 8.196225e-04, 1.318868e-03, 2.908837e-03, 5.669037e-03,
                                                                                      1.285600e-02, 9.128714e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.567056e-03, 2.257787e-03, 3.625568e-03, 8.009342e-03, 1.545416e-02,
                                                                                      3.579453e-02, 2.686405e-01],
                                                                                     [3.567056e-03, 2.257787e-03, 3.625568e-03, 8.009342e-03, 1.545416e-02,
                                                                                      3.579453e-02, 2.686405e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.339417e-03, 2.122130e-03, 3.435606e-03, 7.532806e-03, 1.479704e-02,
                                                                                      3.390926e-02, 2.409465e-01],
                                                                                     [3.339417e-03, 2.122130e-03, 3.435606e-03, 7.532806e-03, 1.479704e-02,
                                                                                      3.390926e-02, 2.409465e-01],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
