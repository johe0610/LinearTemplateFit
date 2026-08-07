
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 6.500000e+01, 1.150000e+02, 1.700000e+02, 2.500000e+02,
                                                                                   4.000000e+02, 7.000000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 6.500000e+01, 1.150000e+02, 1.700000e+02, 2.500000e+02,
                                                                                   4.000000e+02, 7.000000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 6.500000e+01, 1.150000e+02, 1.700000e+02, 2.500000e+02,
                                                                                   4.000000e+02, 7.000000e+02, 9.500000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 9.000000e+01, 1.400000e+02, 2.000000e+02,
                                                                                   3.000000e+02, 5.000000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 9.000000e+01, 1.400000e+02, 2.000000e+02,
                                                                                   3.000000e+02, 5.000000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 9.000000e+01, 1.400000e+02, 2.000000e+02,
                                                                                   3.000000e+02, 5.000000e+02, 9.000000e+02, 1.000000e+03],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [7.187140e-02, 2.469596e-01, 2.246413e-02, 1.729008e-03, 3.012300e-04,
                                                                                   3.151453e-05, 1.238787e-06, 9.034353e-08],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [7.120475e-02, 2.484746e-01, 2.025550e-02, 8.457075e-04, 8.588251e-05,
                                                                                   4.820881e-06, 7.535355e-08, 0.000000e+00],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.942174e-02, 2.390177e-01, 2.006467e-02, 8.420718e-04, 8.717824e-05,
                                                                                   4.679954e-06, 1.239375e-07, 0.000000e+00],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.205943e-05, 5.303502e-05, 1.599916e-05, 4.064195e-06, 1.319264e-06,
                                                                                      3.034788e-07, 4.385174e-08, 2.258589e-08],
                                                                                     [3.205943e-05, 5.303502e-05, 1.599916e-05, 4.064195e-06, 1.319264e-06,
                                                                                      3.034788e-07, 4.385174e-08, 2.258589e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.233776e-05, 1.038930e-04, 2.967604e-05, 5.580544e-06, 1.384950e-06,
                                                                                      2.357621e-07, 2.013909e-08, 0.000000e+00],
                                                                                     [6.233776e-05, 1.038930e-04, 2.967604e-05, 5.580544e-06, 1.384950e-06,
                                                                                      2.357621e-07, 2.013909e-08, 0.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.906830e-05, 9.779263e-05, 2.834603e-05, 5.326748e-06, 1.341512e-06,
                                                                                      2.243534e-07, 2.478751e-08, 0.000000e+00],
                                                                                     [5.906830e-05, 9.779263e-05, 2.834603e-05, 5.326748e-06, 1.341512e-06,
                                                                                      2.243534e-07, 2.478751e-08, 0.000000e+00],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.907244e-01, 1.006135e+00, 9.016819e-01, 4.891287e-01, 2.851061e-01,
                                                                                   1.529733e-01, 6.082850e-02, 0.000000e+00],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.659161e-01, 9.678413e-01, 8.931870e-01, 4.870260e-01, 2.894076e-01,
                                                                                   1.485015e-01, 1.000475e-01, 0.000000e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.460666e-04, 2.147518e-04, 7.122092e-04, 2.350594e-03, 4.379590e-03,
                                                                                      9.629806e-03, 3.539893e-02, 2.500001e-01],
                                                                                     [4.460666e-04, 2.147518e-04, 7.122092e-04, 2.350594e-03, 4.379590e-03,
                                                                                      9.629806e-03, 3.539893e-02, 2.500001e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.673514e-04, 4.206882e-04, 1.321041e-03, 3.227599e-03, 4.597650e-03,
                                                                                      7.481060e-03, 1.625710e-02, 0.000000e+00],
                                                                                     [8.673514e-04, 4.206882e-04, 1.321041e-03, 3.227599e-03, 4.597650e-03,
                                                                                      7.481060e-03, 1.625710e-02, 0.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.218610e-04, 3.959863e-04, 1.261835e-03, 3.080812e-03, 4.453448e-03,
                                                                                      7.119046e-03, 2.000950e-02, 0.000000e+00],
                                                                                     [8.218610e-04, 3.959863e-04, 1.261835e-03, 3.080812e-03, 4.453448e-03,
                                                                                      7.119046e-03, 2.000950e-02, 0.000000e+00],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
