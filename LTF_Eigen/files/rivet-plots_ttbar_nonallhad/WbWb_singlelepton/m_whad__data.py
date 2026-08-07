
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.800000e+01, 7.600000e+01, 8.400000e+01, 9.200000e+01],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.800000e+01, 7.600000e+01, 8.400000e+01, 9.200000e+01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.800000e+01, 7.600000e+01, 8.400000e+01, 9.200000e+01],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.400000e+01, 7.200000e+01, 8.000000e+01, 8.800000e+01, 9.600000e+01],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.400000e+01, 7.200000e+01, 8.000000e+01, 8.800000e+01, 9.600000e+01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.400000e+01, 7.200000e+01, 8.000000e+01, 8.800000e+01, 9.600000e+01],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.775812e-01, 7.300767e-01, 6.957009e-01, 1.925065e-01],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.719637e-01, 7.315212e-01, 6.946913e-01, 1.868617e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.609064e-01, 7.071996e-01, 6.752668e-01, 1.793527e-01],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                     [4.000000e+00, 4.000000e+00, 4.000000e+00, 4.000000e+00],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.405100e-04, 2.279643e-04, 2.228430e-04, 1.171408e-04],
                                                                                     [1.405100e-04, 2.279643e-04, 2.228430e-04, 1.171408e-04],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.716124e-04, 4.456645e-04, 4.349290e-04, 2.254300e-04],
                                                                                     [2.716124e-04, 4.456645e-04, 4.349290e-04, 2.254300e-04],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.553161e-04, 4.205310e-04, 4.115181e-04, 2.119871e-04],
                                                                                     [2.553161e-04, 4.205310e-04, 4.115181e-04, 2.119871e-04],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.797627e-01, 1.001979e+00, 9.985488e-01, 9.706774e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.399282e-01, 9.686648e-01, 9.706280e-01, 9.316709e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.061942e-04, 3.122471e-04, 3.203144e-04, 6.085031e-04],
                                                                                     [5.061942e-04, 3.122471e-04, 3.203144e-04, 6.085031e-04],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.784971e-04, 6.104352e-04, 6.251666e-04, 1.171025e-03],
                                                                                     [9.784971e-04, 6.104352e-04, 6.251666e-04, 1.171025e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.197889e-04, 5.760093e-04, 5.915158e-04, 1.101195e-03],
                                                                                     [9.197889e-04, 5.760093e-04, 5.915158e-04, 1.101195e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
