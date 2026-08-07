
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 3.800000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 3.800000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 3.800000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   6.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   6.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.410183e-02, 5.692952e-03, 1.530798e-03, 6.556533e-04, 6.443131e-05],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.295017e-02, 5.531891e-03, 1.468960e-03, 6.358853e-04, 6.169207e-05],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.127078e-02, 5.386617e-03, 1.443034e-03, 6.238868e-04, 5.951247e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.492188e-05, 1.042596e-05, 5.002067e-06, 3.056747e-06, 2.886561e-07],
                                                                                     [3.492188e-05, 1.042596e-05, 5.002067e-06, 3.056747e-06, 2.886561e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.759146e-05, 2.007615e-05, 9.563809e-06, 5.884137e-06, 5.516179e-07],
                                                                                     [6.759146e-05, 2.007615e-05, 9.563809e-06, 5.884137e-06, 5.516179e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.399443e-05, 1.902248e-05, 9.104575e-06, 5.589769e-06, 5.199078e-07],
                                                                                     [6.399443e-05, 1.902248e-05, 9.104575e-06, 5.589769e-06, 5.199078e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.820339e-01, 9.717087e-01, 9.596041e-01, 9.698499e-01, 9.574859e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.558351e-01, 9.461905e-01, 9.426678e-01, 9.515499e-01, 9.236576e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.447876e-04, 1.831380e-03, 3.267621e-03, 4.662139e-03, 4.480059e-03],
                                                                                     [5.447876e-04, 1.831380e-03, 3.267621e-03, 4.662139e-03, 4.480059e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.054439e-03, 3.526492e-03, 6.247597e-03, 8.974464e-03, 8.561333e-03],
                                                                                     [1.054439e-03, 3.526492e-03, 6.247597e-03, 8.974464e-03, 8.561333e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.983245e-04, 3.341409e-03, 5.947601e-03, 8.525495e-03, 8.069179e-03],
                                                                                     [9.983245e-04, 3.341409e-03, 5.947601e-03, 8.525495e-03, 8.069179e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
