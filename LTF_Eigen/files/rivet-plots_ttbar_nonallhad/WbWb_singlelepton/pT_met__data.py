
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 8.500000e+01, 1.650000e+02, 2.650000e+02, 4.100000e+02,
                                                                                   7.000000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 8.500000e+01, 1.650000e+02, 2.650000e+02, 4.100000e+02,
                                                                                   7.000000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 8.500000e+01, 1.650000e+02, 2.650000e+02, 4.100000e+02,
                                                                                   7.000000e+02, 9.500000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.200000e+02, 2.100000e+02, 3.200000e+02,
                                                                                   5.000000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.200000e+02, 2.100000e+02, 3.200000e+02,
                                                                                   5.000000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.200000e+02, 2.100000e+02, 3.200000e+02,
                                                                                   5.000000e+02, 9.000000e+02, 1.000000e+03],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.453636e-01, 1.038068e-01, 1.814491e-02, 2.422782e-03, 2.638656e-04,
                                                                                   1.232575e-05, 4.797318e-07],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.467526e-01, 1.021347e-01, 1.734735e-02, 2.270672e-03, 2.375740e-04,
                                                                                   1.019726e-05, 4.089422e-07],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.402881e-01, 9.915473e-02, 1.712670e-02, 2.245484e-03, 2.407194e-04,
                                                                                   1.041134e-05, 4.164500e-07],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 3.500000e+01, 4.500000e+01, 5.500000e+01, 9.000000e+01,
                                                                                      2.000000e+02, 5.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.064838e-05, 2.906409e-05, 1.076838e-05, 3.594002e-06, 9.390924e-07,
                                                                                      1.376641e-07, 5.558734e-08],
                                                                                     [4.064838e-05, 2.906409e-05, 1.076838e-05, 3.594002e-06, 9.390924e-07,
                                                                                      1.376641e-07, 5.558734e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.975986e-05, 5.631044e-05, 2.057201e-05, 6.797763e-06, 1.746922e-06,
                                                                                      2.463017e-07, 9.381779e-08],
                                                                                     [7.975986e-05, 5.631044e-05, 2.057201e-05, 6.797763e-06, 1.746922e-06,
                                                                                      2.463017e-07, 9.381779e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.484738e-05, 5.324432e-05, 1.960878e-05, 6.488825e-06, 1.686313e-06,
                                                                                      2.389025e-07, 9.087685e-08],
                                                                                     [7.484738e-05, 5.324432e-05, 1.960878e-05, 6.488825e-06, 1.686313e-06,
                                                                                      2.389025e-07, 9.087685e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.009555e+00, 9.838922e-01, 9.560450e-01, 9.372168e-01, 9.003599e-01,
                                                                                   8.273136e-01, 8.524392e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.650841e-01, 9.551853e-01, 9.438845e-01, 9.268205e-01, 9.122803e-01,
                                                                                   8.446821e-01, 8.680892e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.796325e-04, 2.799825e-04, 5.934656e-04, 1.483419e-03, 3.558980e-03,
                                                                                      1.116882e-02, 1.158717e-01],
                                                                                     [2.796325e-04, 2.799825e-04, 5.934656e-04, 1.483419e-03, 3.558980e-03,
                                                                                      1.116882e-02, 1.158717e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.486921e-04, 5.424543e-04, 1.133762e-03, 2.805768e-03, 6.620499e-03,
                                                                                      1.998269e-02, 1.955630e-01],
                                                                                     [5.486921e-04, 5.424543e-04, 1.133762e-03, 2.805768e-03, 6.620499e-03,
                                                                                      1.998269e-02, 1.955630e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.148977e-04, 5.129175e-04, 1.080677e-03, 2.678254e-03, 6.390803e-03,
                                                                                      1.938239e-02, 1.894326e-01],
                                                                                     [5.148977e-04, 5.129175e-04, 1.080677e-03, 2.678254e-03, 6.390803e-03,
                                                                                      1.938239e-02, 1.894326e-01],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
