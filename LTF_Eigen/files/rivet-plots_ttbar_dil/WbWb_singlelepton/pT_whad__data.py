
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.600000e+02,
                                                                                   3.000000e+02, 3.450000e+02, 4.000000e+02, 4.750000e+02, 6.000000e+02,
                                                                                   7.900000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.600000e+02,
                                                                                   3.000000e+02, 3.450000e+02, 4.000000e+02, 4.750000e+02, 6.000000e+02,
                                                                                   7.900000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.600000e+02,
                                                                                   3.000000e+02, 3.450000e+02, 4.000000e+02, 4.750000e+02, 6.000000e+02,
                                                                                   7.900000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   2.800000e+02, 3.200000e+02, 3.700000e+02, 4.300000e+02, 5.200000e+02,
                                                                                   6.800000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   2.800000e+02, 3.200000e+02, 3.700000e+02, 4.300000e+02, 5.200000e+02,
                                                                                   6.800000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   2.800000e+02, 3.200000e+02, 3.700000e+02, 4.300000e+02, 5.200000e+02,
                                                                                   6.800000e+02, 9.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.633397e-03, 2.117976e-03, 1.103603e-03, 6.381002e-04, 3.879295e-04,
                                                                                   2.534092e-04, 1.742261e-04, 1.191891e-04, 7.539175e-05, 3.675484e-05,
                                                                                   1.176452e-05],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.726343e-03, 2.174998e-03, 1.136723e-03, 6.561587e-04, 3.959230e-04,
                                                                                   2.619171e-04, 1.817979e-04, 1.222038e-04, 7.824586e-05, 3.750129e-05,
                                                                                   1.232104e-05],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.561958e-03, 2.085766e-03, 1.094412e-03, 6.226875e-04, 3.776613e-04,
                                                                                   2.502252e-04, 1.717251e-04, 1.183382e-04, 7.639227e-05, 3.611177e-05,
                                                                                   1.158950e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.044464e-06, 2.735827e-06, 1.974905e-06, 1.501339e-06, 1.170359e-06,
                                                                                      9.449319e-07, 7.009203e-07, 5.288761e-07, 3.433143e-07, 1.796098e-07,
                                                                                      8.662322e-08],
                                                                                     [4.044464e-06, 2.735827e-06, 1.974905e-06, 1.501339e-06, 1.170359e-06,
                                                                                      9.449319e-07, 7.009203e-07, 5.288761e-07, 3.433143e-07, 1.796098e-07,
                                                                                      8.662322e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.110866e-05, 7.542338e-06, 5.452466e-06, 4.139671e-06, 3.212557e-06,
                                                                                      2.613240e-06, 1.946585e-06, 1.457153e-06, 9.514511e-07, 4.932233e-07,
                                                                                      2.410334e-07],
                                                                                     [1.110866e-05, 7.542338e-06, 5.452466e-06, 4.139671e-06, 3.212557e-06,
                                                                                      2.613240e-06, 1.946585e-06, 1.457153e-06, 9.514511e-07, 4.932233e-07,
                                                                                      2.410334e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.047368e-05, 7.087403e-06, 5.128291e-06, 3.873192e-06, 3.010280e-06,
                                                                                      2.450670e-06, 1.815967e-06, 1.376462e-06, 9.017259e-07, 4.648442e-07,
                                                                                      2.242897e-07],
                                                                                     [1.047368e-05, 7.087403e-06, 5.128291e-06, 3.873192e-06, 3.010280e-06,
                                                                                      2.450670e-06, 1.815967e-06, 1.376462e-06, 9.017259e-07, 4.648442e-07,
                                                                                      2.242897e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.020060e+00, 1.026923e+00, 1.030011e+00, 1.028300e+00, 1.020606e+00,
                                                                                   1.033574e+00, 1.043460e+00, 1.025293e+00, 1.037857e+00, 1.020309e+00,
                                                                                   1.047305e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.845817e-01, 9.847921e-01, 9.916718e-01, 9.758460e-01, 9.735308e-01,
                                                                                   9.874353e-01, 9.856451e-01, 9.928609e-01, 1.013271e+00, 9.825038e-01,
                                                                                   9.851231e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.728939e-04, 1.291718e-03, 1.789507e-03, 2.352826e-03, 3.016937e-03,
                                                                                      3.728878e-03, 4.023050e-03, 4.437286e-03, 4.553738e-03, 4.886698e-03,
                                                                                      7.363090e-03],
                                                                                     [8.728939e-04, 1.291718e-03, 1.789507e-03, 2.352826e-03, 3.016937e-03,
                                                                                      3.728878e-03, 4.023050e-03, 4.437286e-03, 4.553738e-03, 4.886698e-03,
                                                                                      7.363090e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.397520e-03, 3.561106e-03, 4.940605e-03, 6.487494e-03, 8.281291e-03,
                                                                                      1.031233e-02, 1.117275e-02, 1.222556e-02, 1.262010e-02, 1.341927e-02,
                                                                                      2.048816e-02],
                                                                                     [2.397520e-03, 3.561106e-03, 4.940605e-03, 6.487494e-03, 8.281291e-03,
                                                                                      1.031233e-02, 1.117275e-02, 1.222556e-02, 1.262010e-02, 1.341927e-02,
                                                                                      2.048816e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.260475e-03, 3.346309e-03, 4.646862e-03, 6.069881e-03, 7.759864e-03,
                                                                                      9.670801e-03, 1.042305e-02, 1.154856e-02, 1.196054e-02, 1.264716e-02,
                                                                                      1.906493e-02],
                                                                                     [2.260475e-03, 3.346309e-03, 4.646862e-03, 6.069881e-03, 7.759864e-03,
                                                                                      9.670801e-03, 1.042305e-02, 1.154856e-02, 1.196054e-02, 1.264716e-02,
                                                                                      1.906493e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
