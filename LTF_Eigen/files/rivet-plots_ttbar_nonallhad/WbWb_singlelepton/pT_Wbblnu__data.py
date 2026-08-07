
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.250000e+02, 1.850000e+02, 2.800000e+02,
                                                                                   4.100000e+02, 5.800000e+02, 8.300000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.250000e+02, 1.850000e+02, 2.800000e+02,
                                                                                   4.100000e+02, 5.800000e+02, 8.300000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.250000e+02, 1.850000e+02, 2.800000e+02,
                                                                                   4.100000e+02, 5.800000e+02, 8.300000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.500000e+02, 2.200000e+02,
                                                                                   3.400000e+02, 4.800000e+02, 6.800000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.500000e+02, 2.200000e+02,
                                                                                   3.400000e+02, 4.800000e+02, 6.800000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.500000e+02, 2.200000e+02,
                                                                                   3.400000e+02, 4.800000e+02, 6.800000e+02, 9.800000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.093934e-01, 8.814600e-02, 2.066489e-02, 5.192416e-03, 1.149230e-03,
                                                                                   2.943517e-04, 1.061647e-04, 3.412946e-05],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.084793e-01, 8.695231e-02, 2.026596e-02, 5.043417e-03, 1.120947e-03,
                                                                                   2.817124e-04, 1.044150e-04, 3.300337e-05],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.006504e-01, 8.446604e-02, 1.973660e-02, 4.954907e-03, 1.086572e-03,
                                                                                   2.761510e-04, 1.022500e-04, 3.299243e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 2.500000e+01, 3.500000e+01, 6.000000e+01,
                                                                                      7.000000e+01, 1.000000e+02, 1.500000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.882857e-05, 3.171988e-05, 1.537809e-05, 6.518893e-06, 2.340912e-06,
                                                                                      1.094883e-06, 5.496365e-07, 2.539885e-07],
                                                                                     [4.882857e-05, 3.171988e-05, 1.537809e-05, 6.518893e-06, 2.340912e-06,
                                                                                      1.094883e-06, 5.496365e-07, 2.539885e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.515888e-05, 6.152659e-05, 2.974939e-05, 1.254883e-05, 4.515430e-06,
                                                                                      2.092612e-06, 1.063664e-06, 4.877716e-07],
                                                                                     [9.515888e-05, 6.152659e-05, 2.974939e-05, 1.254883e-05, 4.515430e-06,
                                                                                      2.092612e-06, 1.063664e-06, 4.877716e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.959129e-05, 5.819732e-05, 2.817458e-05, 1.194073e-05, 4.267983e-06,
                                                                                      1.989444e-06, 1.010227e-06, 4.678457e-07],
                                                                                     [8.959129e-05, 5.819732e-05, 2.817458e-05, 1.194073e-05, 4.267983e-06,
                                                                                      1.989444e-06, 1.010227e-06, 4.678457e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.956345e-01, 9.864578e-01, 9.806953e-01, 9.713045e-01, 9.753896e-01,
                                                                                   9.570606e-01, 9.835190e-01, 9.670053e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.582461e-01, 9.582515e-01, 9.550789e-01, 9.542585e-01, 9.454783e-01,
                                                                                   9.381668e-01, 9.631262e-01, 9.666848e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.331906e-04, 3.598561e-04, 7.441651e-04, 1.255464e-03, 2.036940e-03,
                                                                                      3.719642e-03, 5.177206e-03, 7.441914e-03],
                                                                                     [2.331906e-04, 3.598561e-04, 7.441651e-04, 1.255464e-03, 2.036940e-03,
                                                                                      3.719642e-03, 5.177206e-03, 7.441914e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.544502e-04, 6.980077e-04, 1.439610e-03, 2.416761e-03, 3.929092e-03,
                                                                                      7.109223e-03, 1.001900e-02, 1.429181e-02],
                                                                                     [4.544502e-04, 6.980077e-04, 1.439610e-03, 2.416761e-03, 3.929092e-03,
                                                                                      7.109223e-03, 1.001900e-02, 1.429181e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.278611e-04, 6.602378e-04, 1.363403e-03, 2.299648e-03, 3.713776e-03,
                                                                                      6.758731e-03, 9.515658e-03, 1.370797e-02],
                                                                                     [4.278611e-04, 6.602378e-04, 1.363403e-03, 2.299648e-03, 3.713776e-03,
                                                                                      6.758731e-03, 9.515658e-03, 1.370797e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
