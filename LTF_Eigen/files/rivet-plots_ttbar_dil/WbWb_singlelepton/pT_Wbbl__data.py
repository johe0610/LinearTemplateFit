
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 9.000000e+01, 1.500000e+02, 2.100000e+02, 2.700000e+02,
                                                                                   3.350000e+02, 4.100000e+02, 5.000000e+02, 6.250000e+02, 8.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 9.000000e+01, 1.500000e+02, 2.100000e+02, 2.700000e+02,
                                                                                   3.350000e+02, 4.100000e+02, 5.000000e+02, 6.250000e+02, 8.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 9.000000e+01, 1.500000e+02, 2.100000e+02, 2.700000e+02,
                                                                                   3.350000e+02, 4.100000e+02, 5.000000e+02, 6.250000e+02, 8.000000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 6.000000e+01, 1.200000e+02, 1.800000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.700000e+02, 4.500000e+02, 5.500000e+02, 7.000000e+02,
                                                                                   9.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 6.000000e+01, 1.200000e+02, 1.800000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.700000e+02, 4.500000e+02, 5.500000e+02, 7.000000e+02,
                                                                                   9.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 6.000000e+01, 1.200000e+02, 1.800000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.700000e+02, 4.500000e+02, 5.500000e+02, 7.000000e+02,
                                                                                   9.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.774010e-03, 2.760592e-03, 1.301989e-03, 4.464878e-04, 1.661807e-04,
                                                                                   7.268702e-05, 3.445852e-05, 1.503947e-05, 5.525671e-06, 1.432519e-06],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.822982e-03, 2.828841e-03, 1.330457e-03, 4.587020e-04, 1.712293e-04,
                                                                                   7.370127e-05, 3.269176e-05, 1.514324e-05, 5.530886e-06, 1.569063e-06],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.732747e-03, 2.710789e-03, 1.293417e-03, 4.493038e-04, 1.652731e-04,
                                                                                   7.342034e-05, 3.326657e-05, 1.464751e-05, 5.738112e-06, 1.512483e-06],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.042178e-06, 2.548322e-06, 1.751286e-06, 1.026914e-06, 6.265649e-07,
                                                                                      3.835946e-07, 2.467764e-07, 1.458220e-07, 7.211310e-08, 3.183975e-08],
                                                                                     [2.042178e-06, 2.548322e-06, 1.751286e-06, 1.026914e-06, 6.265649e-07,
                                                                                      3.835946e-07, 2.467764e-07, 1.458220e-07, 7.211310e-08, 3.183975e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.630258e-06, 7.015635e-06, 4.815256e-06, 2.828015e-06, 1.732568e-06,
                                                                                      1.051175e-06, 6.564913e-07, 3.968564e-07, 1.961554e-07, 9.014035e-08],
                                                                                     [5.630258e-06, 7.015635e-06, 4.815256e-06, 2.828015e-06, 1.732568e-06,
                                                                                      1.051175e-06, 6.564913e-07, 3.968564e-07, 1.961554e-07, 9.014035e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.267897e-06, 6.590746e-06, 4.554329e-06, 2.687207e-06, 1.629665e-06,
                                                                                      1.006535e-06, 6.332979e-07, 3.770153e-07, 1.914814e-07, 8.521727e-08],
                                                                                     [5.267897e-06, 6.590746e-06, 4.554329e-06, 2.687207e-06, 1.629665e-06,
                                                                                      1.006535e-06, 6.332979e-07, 3.770153e-07, 1.914814e-07, 8.521727e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.027605e+00, 1.024723e+00, 1.021865e+00, 1.027356e+00, 1.030380e+00,
                                                                                   1.013954e+00, 9.487279e-01, 1.006900e+00, 1.000944e+00, 1.095317e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.767403e-01, 9.819593e-01, 9.934162e-01, 1.006307e+00, 9.945385e-01,
                                                                                   1.010089e+00, 9.654091e-01, 9.739379e-01, 1.038446e+00, 1.055821e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.151165e-03, 9.231071e-04, 1.345085e-03, 2.299982e-03, 3.770383e-03,
                                                                                      5.277347e-03, 7.161550e-03, 9.695953e-03, 1.305056e-02, 2.222641e-02],
                                                                                     [1.151165e-03, 9.231071e-04, 1.345085e-03, 2.299982e-03, 3.770383e-03,
                                                                                      5.277347e-03, 7.161550e-03, 9.695953e-03, 1.305056e-02, 2.222641e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.173746e-03, 2.541352e-03, 3.698385e-03, 6.333913e-03, 1.042581e-02,
                                                                                      1.446166e-02, 1.905164e-02, 2.638766e-02, 3.549893e-02, 6.292437e-02],
                                                                                     [3.173746e-03, 2.541352e-03, 3.698385e-03, 6.333913e-03, 1.042581e-02,
                                                                                      1.446166e-02, 1.905164e-02, 2.638766e-02, 3.549893e-02, 6.292437e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.969486e-03, 2.387439e-03, 3.497978e-03, 6.018545e-03, 9.806584e-03,
                                                                                      1.384752e-02, 1.837856e-02, 2.506839e-02, 3.465306e-02, 5.948771e-02],
                                                                                     [2.969486e-03, 2.387439e-03, 3.497978e-03, 6.018545e-03, 9.806584e-03,
                                                                                      1.384752e-02, 1.837856e-02, 2.506839e-02, 3.465306e-02, 5.948771e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
