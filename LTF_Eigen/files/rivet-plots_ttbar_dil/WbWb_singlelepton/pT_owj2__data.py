
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 3.800000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 3.800000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 3.800000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   6.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   6.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.900163e-03, 2.300637e-04, 7.499005e-05, 3.140913e-05, 2.710916e-06],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.953281e-03, 2.352236e-04, 7.262280e-05, 3.215993e-05, 2.756523e-06],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.882383e-03, 2.279052e-04, 7.369375e-05, 3.086980e-05, 2.721675e-06],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.200000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.990279e-06, 1.039556e-06, 5.492113e-07, 3.327415e-07, 2.942462e-08],
                                                                                     [2.990279e-06, 1.039556e-06, 5.492113e-07, 3.327415e-07, 2.942462e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.246108e-06, 2.858969e-06, 1.472563e-06, 9.132816e-07, 8.062200e-08],
                                                                                     [8.246108e-06, 2.858969e-06, 1.472563e-06, 9.132816e-07, 8.062200e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.764612e-06, 2.696570e-06, 1.420593e-06, 8.621277e-07, 7.694923e-08],
                                                                                     [7.764612e-06, 2.696570e-06, 1.420593e-06, 8.621277e-07, 7.694923e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.027954e+00, 1.022428e+00, 9.684325e-01, 1.023904e+00, 1.016823e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.906429e-01, 9.906178e-01, 9.827137e-01, 9.828289e-01, 1.003969e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.573696e-03, 4.518557e-03, 7.323789e-03, 1.059378e-02, 1.085412e-02],
                                                                                     [1.573696e-03, 4.518557e-03, 7.323789e-03, 1.059378e-02, 1.085412e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.339685e-03, 1.242686e-02, 1.963678e-02, 2.907695e-02, 2.973976e-02],
                                                                                     [4.339685e-03, 1.242686e-02, 1.963678e-02, 2.907695e-02, 2.973976e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.086287e-03, 1.172097e-02, 1.894375e-02, 2.744832e-02, 2.838496e-02],
                                                                                     [4.086287e-03, 1.172097e-02, 1.894375e-02, 2.744832e-02, 2.838496e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
