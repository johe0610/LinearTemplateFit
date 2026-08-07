
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.100000e+02, 4.800000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.100000e+02, 4.800000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.100000e+02, 4.800000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   2.000000e+02, 2.600000e+02, 3.600000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   2.000000e+02, 2.600000e+02, 3.600000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   2.000000e+02, 2.600000e+02, 3.600000e+02, 6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.744240e-03, 7.004933e-04, 1.875070e-04, 8.820980e-05, 5.474629e-05,
                                                                                   3.413108e-05, 1.642287e-05, 4.524211e-06],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.859805e-03, 7.223119e-04, 1.948147e-04, 8.905068e-05, 5.375608e-05,
                                                                                   3.495708e-05, 1.681044e-05, 4.716939e-06],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.680492e-03, 6.950335e-04, 1.862680e-04, 8.671512e-05, 5.477408e-05,
                                                                                   3.392246e-05, 1.596383e-05, 4.278089e-06],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.725665e-06, 1.815862e-06, 8.690564e-07, 5.572226e-07, 4.387948e-07,
                                                                                      2.829137e-07, 1.520872e-07, 5.145774e-08],
                                                                                     [4.725665e-06, 1.815862e-06, 8.690564e-07, 5.572226e-07, 4.387948e-07,
                                                                                      2.829137e-07, 1.520872e-07, 5.145774e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.300836e-05, 5.013834e-06, 2.409111e-06, 1.522931e-06, 1.182086e-06,
                                                                                      7.787405e-07, 4.185551e-07, 1.433275e-07],
                                                                                     [1.300836e-05, 5.013834e-06, 2.409111e-06, 1.522931e-06, 1.182086e-06,
                                                                                      7.787405e-07, 4.185551e-07, 1.433275e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.224932e-05, 4.720123e-06, 2.260122e-06, 1.441834e-06, 1.146096e-06,
                                                                                      7.358754e-07, 3.907581e-07, 1.305410e-07],
                                                                                     [1.224932e-05, 4.720123e-06, 2.260122e-06, 1.441834e-06, 1.146096e-06,
                                                                                      7.358754e-07, 3.907581e-07, 1.305410e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.024359e+00, 1.031147e+00, 1.038973e+00, 1.009533e+00, 9.819127e-01,
                                                                                   1.024201e+00, 1.023599e+00, 1.042599e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.865631e-01, 9.922058e-01, 9.933922e-01, 9.830554e-01, 1.000508e+00,
                                                                                   9.938877e-01, 9.720487e-01, 9.455989e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.960847e-04, 2.592262e-03, 4.634794e-03, 6.317015e-03, 8.015060e-03,
                                                                                      8.289035e-03, 9.260696e-03, 1.137386e-02],
                                                                                     [9.960847e-04, 2.592262e-03, 4.634794e-03, 6.317015e-03, 8.015060e-03,
                                                                                      8.289035e-03, 9.260696e-03, 1.137386e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.741927e-03, 7.157576e-03, 1.284811e-02, 1.726487e-02, 2.159208e-02,
                                                                                      2.281617e-02, 2.548611e-02, 3.168011e-02],
                                                                                     [2.741927e-03, 7.157576e-03, 1.284811e-02, 1.726487e-02, 2.159208e-02,
                                                                                      2.281617e-02, 2.548611e-02, 3.168011e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.581935e-03, 6.738284e-03, 1.205353e-02, 1.634551e-02, 2.093468e-02,
                                                                                      2.156027e-02, 2.379353e-02, 2.885387e-02],
                                                                                     [2.581935e-03, 6.738284e-03, 1.205353e-02, 1.634551e-02, 2.093468e-02,
                                                                                      2.156027e-02, 2.379353e-02, 2.885387e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
