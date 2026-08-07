
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [7.694628e-02, 2.420012e-01, 2.667213e-01, 2.021325e-01, 8.945254e-02,
                                                                                   2.133767e-02, 2.630329e-03],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.010367e-02, 2.517798e-01, 2.762757e-01, 2.063639e-01, 9.038295e-02,
                                                                                   2.122119e-02, 2.566057e-03],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [7.470892e-02, 2.348662e-01, 2.606052e-01, 1.986831e-01, 8.945942e-02,
                                                                                   2.181784e-02, 2.669059e-03],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.644834e-04, 3.369630e-04, 3.538999e-04, 2.388270e-04, 1.186226e-04,
                                                                                      5.035710e-05, 1.259196e-05],
                                                                                     [1.644834e-04, 3.369630e-04, 3.538999e-04, 2.388270e-04, 1.186226e-04,
                                                                                      5.035710e-05, 1.259196e-05],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.562748e-04, 9.346500e-04, 9.794450e-04, 6.564018e-04, 3.242451e-04,
                                                                                      1.366656e-04, 3.399217e-05],
                                                                                     [4.562748e-04, 9.346500e-04, 9.794450e-04, 6.564018e-04, 3.242451e-04,
                                                                                      1.366656e-04, 3.399217e-05],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.229326e-04, 8.664656e-04, 9.129313e-04, 6.178728e-04, 3.095179e-04,
                                                                                      1.328124e-04, 3.324377e-05],
                                                                                     [4.229326e-04, 8.664656e-04, 9.129313e-04, 6.178728e-04, 3.095179e-04,
                                                                                      1.328124e-04, 3.324377e-05],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.041034e+00, 1.040407e+00, 1.035822e+00, 1.020934e+00, 1.010401e+00,
                                                                                   9.945411e-01, 9.755650e-01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.709231e-01, 9.705167e-01, 9.770693e-01, 9.829350e-01, 1.000077e+00,
                                                                                   1.022503e+00, 1.014724e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.137639e-03, 1.392402e-03, 1.326853e-03, 1.181537e-03, 1.326095e-03,
                                                                                      2.360009e-03, 4.787219e-03],
                                                                                     [2.137639e-03, 1.392402e-03, 1.326853e-03, 1.181537e-03, 1.326095e-03,
                                                                                      2.360009e-03, 4.787219e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.929784e-03, 3.862171e-03, 3.672166e-03, 3.247384e-03, 3.624772e-03,
                                                                                      6.404898e-03, 1.292316e-02],
                                                                                     [5.929784e-03, 3.862171e-03, 3.672166e-03, 3.247384e-03, 3.624772e-03,
                                                                                      6.404898e-03, 1.292316e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.496466e-03, 3.580419e-03, 3.422791e-03, 3.056771e-03, 3.460135e-03,
                                                                                      6.224316e-03, 1.263864e-02],
                                                                                     [5.496466e-03, 3.580419e-03, 3.422791e-03, 3.056771e-03, 3.460135e-03,
                                                                                      6.224316e-03, 1.263864e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
