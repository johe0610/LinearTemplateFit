
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.000000e+01, 9.000000e+01, 1.300000e+02, 1.700000e+02, 2.100000e+02,
                                                                                   2.600000e+02, 3.200000e+02, 4.000000e+02, 5.250000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.000000e+01, 9.000000e+01, 1.300000e+02, 1.700000e+02, 2.100000e+02,
                                                                                   2.600000e+02, 3.200000e+02, 4.000000e+02, 5.250000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.000000e+01, 9.000000e+01, 1.300000e+02, 1.700000e+02, 2.100000e+02,
                                                                                   2.600000e+02, 3.200000e+02, 4.000000e+02, 5.250000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 7.000000e+01, 1.100000e+02, 1.500000e+02, 1.900000e+02,
                                                                                   2.300000e+02, 2.900000e+02, 3.500000e+02, 4.500000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 7.000000e+01, 1.100000e+02, 1.500000e+02, 1.900000e+02,
                                                                                   2.300000e+02, 2.900000e+02, 3.500000e+02, 4.500000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 7.000000e+01, 1.100000e+02, 1.500000e+02, 1.900000e+02,
                                                                                   2.300000e+02, 2.900000e+02, 3.500000e+02, 4.500000e+02, 6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.589801e-01, 9.771099e-02, 3.459542e-02, 1.243469e-02, 4.760209e-03,
                                                                                   1.649394e-03, 5.110365e-04, 1.380808e-04, 2.284615e-05],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.574575e-01, 9.678155e-02, 3.408383e-02, 1.211061e-02, 4.634197e-03,
                                                                                   1.604006e-03, 4.878925e-04, 1.303949e-04, 2.108004e-05],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.478191e-01, 9.368592e-02, 3.322719e-02, 1.195554e-02, 4.563114e-03,
                                                                                   1.585503e-03, 4.863277e-04, 1.299535e-04, 2.172141e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 3.000000e+01, 5.000000e+01, 7.500000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.066718e-05, 3.732883e-05, 2.226764e-05, 1.339321e-05, 8.322145e-06,
                                                                                      4.024182e-06, 2.256797e-06, 9.176988e-07, 3.095825e-07],
                                                                                     [6.066718e-05, 3.732883e-05, 2.226764e-05, 1.339321e-05, 8.322145e-06,
                                                                                      4.024182e-06, 2.256797e-06, 9.176988e-07, 3.095825e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.181348e-04, 7.255649e-05, 4.317962e-05, 2.582198e-05, 1.603451e-05,
                                                                                      7.762256e-06, 4.313892e-06, 1.743346e-06, 5.848428e-07],
                                                                                     [1.181348e-04, 7.255649e-05, 4.317962e-05, 2.582198e-05, 1.603451e-05,
                                                                                      7.762256e-06, 4.313892e-06, 1.743346e-06, 5.848428e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.112333e-04, 6.851011e-05, 4.090766e-05, 2.461787e-05, 1.528189e-05,
                                                                                      7.393444e-06, 4.129784e-06, 1.671347e-06, 5.691055e-07],
                                                                                     [1.112333e-04, 6.851011e-05, 4.090766e-05, 2.461787e-05, 1.528189e-05,
                                                                                      7.393444e-06, 4.129784e-06, 1.671347e-06, 5.691055e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.941208e-01, 9.904879e-01, 9.852122e-01, 9.739374e-01, 9.735281e-01,
                                                                                   9.724820e-01, 9.547116e-01, 9.443377e-01, 9.226955e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.569040e-01, 9.588064e-01, 9.604505e-01, 9.614667e-01, 9.585953e-01,
                                                                                   9.612640e-01, 9.516496e-01, 9.411410e-01, 9.507689e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.342542e-04, 3.820331e-04, 6.436586e-04, 1.077084e-03, 1.748273e-03,
                                                                                      2.439794e-03, 4.416117e-03, 6.646100e-03, 1.355075e-02],
                                                                                     [2.342542e-04, 3.820331e-04, 6.436586e-04, 1.077084e-03, 1.748273e-03,
                                                                                      2.439794e-03, 4.416117e-03, 6.646100e-03, 1.355075e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.561540e-04, 7.425622e-04, 1.248131e-03, 2.076608e-03, 3.368447e-03,
                                                                                      4.706126e-03, 8.441456e-03, 1.262555e-02, 2.559918e-02],
                                                                                     [4.561540e-04, 7.425622e-04, 1.248131e-03, 2.076608e-03, 3.368447e-03,
                                                                                      4.706126e-03, 8.441456e-03, 1.262555e-02, 2.559918e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.295052e-04, 7.011505e-04, 1.182459e-03, 1.979774e-03, 3.210340e-03,
                                                                                      4.482521e-03, 8.081192e-03, 1.210412e-02, 2.491035e-02],
                                                                                     [4.295052e-04, 7.011505e-04, 1.182459e-03, 1.979774e-03, 3.210340e-03,
                                                                                      4.482521e-03, 8.081192e-03, 1.210412e-02, 2.491035e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
