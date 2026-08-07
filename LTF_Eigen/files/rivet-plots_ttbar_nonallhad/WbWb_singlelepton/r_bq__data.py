
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.255103e+00, 8.792109e+00, 1.327770e+01, 1.066955e+01, 3.464680e+00,
                                                                                   6.074238e-01, 6.919876e-02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.239697e+00, 8.826548e+00, 1.333613e+01, 1.060389e+01, 3.364448e+00,
                                                                                   5.816365e-01, 6.570025e-02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.209337e+00, 8.387519e+00, 1.264496e+01, 1.023415e+01, 3.342121e+00,
                                                                                   5.827799e-01, 6.603481e-02],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.342145e-03, 4.089243e-03, 5.019809e-03, 3.484066e-03, 1.481288e-03,
                                                                                      5.387980e-04, 1.294199e-04],
                                                                                     [1.342145e-03, 4.089243e-03, 5.019809e-03, 3.484066e-03, 1.481288e-03,
                                                                                      5.387980e-04, 1.294199e-04],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.605497e-03, 8.001805e-03, 9.826083e-03, 6.783906e-03, 2.850829e-03,
                                                                                      1.029909e-03, 2.461416e-04],
                                                                                     [2.605497e-03, 8.001805e-03, 9.826083e-03, 6.783906e-03, 2.850829e-03,
                                                                                      1.029909e-03, 2.461416e-04],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.470916e-03, 7.486025e-03, 9.183446e-03, 6.395920e-03, 2.726486e-03,
                                                                                      9.889755e-04, 2.367931e-04],
                                                                                     [2.470916e-03, 7.486025e-03, 9.183446e-03, 6.395920e-03, 2.726486e-03,
                                                                                      9.889755e-04, 2.367931e-04],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.877253e-01, 1.003917e+00, 1.004401e+00, 9.938460e-01, 9.710703e-01,
                                                                                   9.575464e-01, 9.494426e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.635361e-01, 9.539826e-01, 9.523457e-01, 9.591923e-01, 9.646262e-01,
                                                                                   9.594288e-01, 9.542774e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.069350e-03, 4.651038e-04, 3.780631e-04, 3.265429e-04, 4.275396e-04,
                                                                                      8.870215e-04, 1.870263e-03],
                                                                                     [1.069350e-03, 4.651038e-04, 3.780631e-04, 3.265429e-04, 4.275396e-04,
                                                                                      8.870215e-04, 1.870263e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.075923e-03, 9.101121e-04, 7.400441e-04, 6.358193e-04, 8.228261e-04,
                                                                                      1.695536e-03, 3.557023e-03],
                                                                                     [2.075923e-03, 9.101121e-04, 7.400441e-04, 6.358193e-04, 8.228261e-04,
                                                                                      1.695536e-03, 3.557023e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.968696e-03, 8.514482e-04, 6.916443e-04, 5.994555e-04, 7.869373e-04,
                                                                                      1.628147e-03, 3.421927e-03],
                                                                                     [1.968696e-03, 8.514482e-04, 6.916443e-04, 5.994555e-04, 7.869373e-04,
                                                                                      1.628147e-03, 3.421927e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
