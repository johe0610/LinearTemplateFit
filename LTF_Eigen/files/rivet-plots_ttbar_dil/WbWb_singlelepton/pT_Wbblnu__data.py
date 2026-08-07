
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.250000e+02, 1.850000e+02, 2.800000e+02,
                                                                                   4.100000e+02, 5.800000e+02, 8.300000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.250000e+02, 1.850000e+02, 2.800000e+02,
                                                                                   4.100000e+02, 5.800000e+02, 8.300000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.250000e+02, 1.850000e+02, 2.800000e+02,
                                                                                   4.100000e+02, 5.800000e+02, 8.300000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.500000e+02, 2.200000e+02,
                                                                                   3.400000e+02, 4.800000e+02, 6.800000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.500000e+02, 2.200000e+02,
                                                                                   3.400000e+02, 4.800000e+02, 6.800000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.500000e+02, 2.200000e+02,
                                                                                   3.400000e+02, 4.800000e+02, 6.800000e+02, 9.800000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.546267e-03, 2.338420e-03, 6.315252e-04, 1.863952e-04, 4.944217e-05,
                                                                                   1.231692e-05, 3.353047e-06, 8.499866e-07],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.654904e-03, 2.400939e-03, 6.465925e-04, 1.912919e-04, 5.047625e-05,
                                                                                   1.220685e-05, 3.360995e-06, 8.389890e-07],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.470628e-03, 2.306271e-03, 6.253286e-04, 1.849319e-04, 4.880821e-05,
                                                                                   1.213234e-05, 3.287376e-06, 8.747013e-07],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.582215e-06, 2.570999e-06, 1.336530e-06, 6.135529e-07, 2.411622e-07,
                                                                                      1.113996e-07, 4.851331e-08, 1.996223e-08],
                                                                                     [3.582215e-06, 2.570999e-06, 1.336530e-06, 6.135529e-07, 2.411622e-07,
                                                                                      1.113996e-07, 4.851331e-08, 1.996223e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.858183e-06, 7.083162e-06, 3.682991e-06, 1.689419e-06, 6.627651e-07,
                                                                                      3.016021e-07, 1.325403e-07, 5.382119e-08],
                                                                                     [9.858183e-06, 7.083162e-06, 3.682991e-06, 1.689419e-06, 6.627651e-07,
                                                                                      3.016021e-07, 1.325403e-07, 5.382119e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.272653e-06, 6.661515e-06, 3.469952e-06, 1.594190e-06, 6.251150e-07,
                                                                                      2.882094e-07, 1.256021e-07, 5.274649e-08],
                                                                                     [9.272653e-06, 6.661515e-06, 3.469952e-06, 1.594190e-06, 6.251150e-07,
                                                                                      2.882094e-07, 1.256021e-07, 5.274649e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.023896e+00, 1.026736e+00, 1.023859e+00, 1.026271e+00, 1.020915e+00,
                                                                                   9.910635e-01, 1.002370e+00, 9.870614e-01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.833624e-01, 9.862518e-01, 9.901879e-01, 9.921495e-01, 9.871777e-01,
                                                                                   9.850141e-01, 9.804145e-01, 1.029077e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.879465e-04, 1.099460e-03, 2.116353e-03, 3.291678e-03, 4.877662e-03,
                                                                                      9.044436e-03, 1.446843e-02, 2.348535e-02],
                                                                                     [7.879465e-04, 1.099460e-03, 2.116353e-03, 3.291678e-03, 4.877662e-03,
                                                                                      9.044436e-03, 1.446843e-02, 2.348535e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.168413e-03, 3.029038e-03, 5.831899e-03, 9.063640e-03, 1.340485e-02,
                                                                                      2.448681e-02, 3.952832e-02, 6.332005e-02],
                                                                                     [2.168413e-03, 3.029038e-03, 5.831899e-03, 9.063640e-03, 1.340485e-02,
                                                                                      2.448681e-02, 3.952832e-02, 6.332005e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.039619e-03, 2.848725e-03, 5.494558e-03, 8.552742e-03, 1.264336e-02,
                                                                                      2.339947e-02, 3.745909e-02, 6.205567e-02],
                                                                                     [2.039619e-03, 2.848725e-03, 5.494558e-03, 8.552742e-03, 1.264336e-02,
                                                                                      2.339947e-02, 3.745909e-02, 6.205567e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
