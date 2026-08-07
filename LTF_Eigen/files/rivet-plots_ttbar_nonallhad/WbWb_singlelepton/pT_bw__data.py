
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 1.200000e+02, 2.000000e+02, 2.800000e+02, 3.800000e+02,
                                                                                   5.200000e+02, 7.500000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 1.200000e+02, 2.000000e+02, 2.800000e+02, 3.800000e+02,
                                                                                   5.200000e+02, 7.500000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 1.200000e+02, 2.000000e+02, 2.800000e+02, 3.800000e+02,
                                                                                   5.200000e+02, 7.500000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 8.000000e+01, 1.600000e+02, 2.400000e+02, 3.200000e+02,
                                                                                   4.400000e+02, 6.000000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 8.000000e+01, 1.600000e+02, 2.400000e+02, 3.200000e+02,
                                                                                   4.400000e+02, 6.000000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 8.000000e+01, 1.600000e+02, 2.400000e+02, 3.200000e+02,
                                                                                   4.400000e+02, 6.000000e+02, 9.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.047490e-02, 8.520512e-02, 5.750646e-02, 2.480166e-02, 8.439578e-03,
                                                                                   2.093965e-03, 3.055684e-04],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.015691e-02, 8.491077e-02, 5.694760e-02, 2.458904e-02, 8.300277e-03,
                                                                                   2.037837e-03, 2.899211e-04],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.960276e-02, 8.113169e-02, 5.527194e-02, 2.393113e-02, 8.160794e-03,
                                                                                   2.016177e-03, 2.898269e-04],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 6.000000e+01,
                                                                                      8.000000e+01, 1.500000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.203193e-05, 2.457074e-05, 2.023630e-05, 1.334097e-05, 6.390932e-06,
                                                                                      2.782460e-06, 7.857338e-07],
                                                                                     [1.203193e-05, 2.457074e-05, 2.023630e-05, 1.334097e-05, 6.390932e-06,
                                                                                      2.782460e-06, 7.857338e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.331363e-05, 4.790403e-05, 3.932887e-05, 2.594587e-05, 1.238301e-05,
                                                                                      5.365882e-06, 1.497256e-06],
                                                                                     [2.331363e-05, 4.790403e-05, 3.932887e-05, 2.594587e-05, 1.238301e-05,
                                                                                      5.365882e-06, 1.497256e-06],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.206464e-05, 4.493989e-05, 3.718542e-05, 2.455670e-05, 1.178080e-05,
                                                                                      5.120281e-06, 1.438761e-06],
                                                                                     [2.206464e-05, 4.493989e-05, 3.718542e-05, 2.455670e-05, 1.178080e-05,
                                                                                      5.120281e-06, 1.438761e-06],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.844693e-01, 9.965454e-01, 9.902818e-01, 9.914272e-01, 9.834943e-01,
                                                                                   9.731953e-01, 9.487928e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.574044e-01, 9.521927e-01, 9.611431e-01, 9.649003e-01, 9.669671e-01,
                                                                                   9.628513e-01, 9.484845e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.876429e-04, 2.883716e-04, 3.518961e-04, 5.379063e-04, 7.572573e-04,
                                                                                      1.328800e-03, 2.571384e-03],
                                                                                     [5.876429e-04, 2.883716e-04, 3.518961e-04, 5.379063e-04, 7.572573e-04,
                                                                                      1.328800e-03, 2.571384e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.138644e-03, 5.622201e-04, 6.839035e-04, 1.046134e-03, 1.467255e-03,
                                                                                      2.562546e-03, 4.899905e-03],
                                                                                     [1.138644e-03, 5.622201e-04, 6.839035e-04, 1.046134e-03, 1.467255e-03,
                                                                                      2.562546e-03, 4.899905e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.077643e-03, 5.274318e-04, 6.466303e-04, 9.901232e-04, 1.395899e-03,
                                                                                      2.445256e-03, 4.708474e-03],
                                                                                     [1.077643e-03, 5.274318e-04, 6.466303e-04, 9.901232e-04, 1.395899e-03,
                                                                                      2.445256e-03, 4.708474e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
