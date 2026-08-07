
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 1.900000e+02, 2.600000e+02,
                                                                                   3.500000e+02, 4.700000e+02, 6.300000e+02, 8.350000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 1.900000e+02, 2.600000e+02,
                                                                                   3.500000e+02, 4.700000e+02, 6.300000e+02, 8.350000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 1.900000e+02, 2.600000e+02,
                                                                                   3.500000e+02, 4.700000e+02, 6.300000e+02, 8.350000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.200000e+02,
                                                                                   3.000000e+02, 4.000000e+02, 5.400000e+02, 7.200000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.200000e+02,
                                                                                   3.000000e+02, 4.000000e+02, 5.400000e+02, 7.200000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.200000e+02,
                                                                                   3.000000e+02, 4.000000e+02, 5.400000e+02, 7.200000e+02, 9.500000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.040704e-04, 1.233311e-03, 2.011825e-03, 1.551406e-03, 7.480429e-04,
                                                                                   2.677257e-04, 8.865866e-05, 2.894686e-05, 8.314776e-06],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.141291e-04, 1.257828e-03, 2.061088e-03, 1.596128e-03, 7.663215e-04,
                                                                                   2.730887e-04, 8.910520e-05, 3.041691e-05, 8.556563e-06],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.993628e-04, 1.207033e-03, 1.974820e-03, 1.532701e-03, 7.334256e-04,
                                                                                   2.699071e-04, 8.900254e-05, 2.852656e-05, 8.256234e-06],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 3.000000e+01, 4.000000e+01,
                                                                                      5.000000e+01, 7.000000e+01, 9.000000e+01, 1.150000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.258122e-07, 1.863900e-06, 2.174185e-06, 1.911152e-06, 1.150713e-06,
                                                                                      6.165146e-07, 2.997732e-07, 1.508207e-07, 7.141489e-08],
                                                                                     [9.258122e-07, 1.863900e-06, 2.174185e-06, 1.911152e-06, 1.150713e-06,
                                                                                      6.165146e-07, 2.997732e-07, 1.508207e-07, 7.141489e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.559845e-06, 5.119553e-06, 5.984508e-06, 5.272032e-06, 3.166945e-06,
                                                                                      1.693969e-06, 8.181760e-07, 4.207284e-07, 1.970221e-07],
                                                                                     [2.559845e-06, 5.119553e-06, 5.984508e-06, 5.272032e-06, 3.166945e-06,
                                                                                      1.693969e-06, 8.181760e-07, 4.207284e-07, 1.970221e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.396634e-06, 4.813369e-06, 5.622433e-06, 4.957334e-06, 2.973003e-06,
                                                                                      1.614135e-06, 7.835081e-07, 3.919776e-07, 1.857285e-07],
                                                                                     [2.396634e-06, 4.813369e-06, 5.622433e-06, 4.957334e-06, 2.973003e-06,
                                                                                      1.614135e-06, 7.835081e-07, 3.919776e-07, 1.857285e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.033080e+00, 1.019879e+00, 1.024487e+00, 1.028827e+00, 1.024435e+00,
                                                                                   1.020032e+00, 1.005037e+00, 1.050784e+00, 1.029079e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.845181e-01, 9.786931e-01, 9.816063e-01, 9.879432e-01, 9.804593e-01,
                                                                                   1.008148e+00, 1.003879e+00, 9.854803e-01, 9.929593e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.044730e-03, 1.511298e-03, 1.080703e-03, 1.231884e-03, 1.538298e-03,
                                                                                      2.302785e-03, 3.381206e-03, 5.210261e-03, 8.588913e-03],
                                                                                     [3.044730e-03, 1.511298e-03, 1.080703e-03, 1.231884e-03, 1.538298e-03,
                                                                                      2.302785e-03, 3.381206e-03, 5.210261e-03, 8.588913e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.418593e-03, 4.151064e-03, 2.974666e-03, 3.398228e-03, 4.233641e-03,
                                                                                      6.327256e-03, 9.228382e-03, 1.453451e-02, 2.369542e-02],
                                                                                     [8.418593e-03, 4.151064e-03, 2.974666e-03, 3.398228e-03, 4.233641e-03,
                                                                                      6.327256e-03, 9.228382e-03, 1.453451e-02, 2.369542e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.881839e-03, 3.902802e-03, 2.794693e-03, 3.195381e-03, 3.974375e-03,
                                                                                      6.029063e-03, 8.837356e-03, 1.354128e-02, 2.233716e-02],
                                                                                     [7.881839e-03, 3.902802e-03, 2.794693e-03, 3.195381e-03, 3.974375e-03,
                                                                                      6.029063e-03, 8.837356e-03, 1.354128e-02, 2.233716e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
