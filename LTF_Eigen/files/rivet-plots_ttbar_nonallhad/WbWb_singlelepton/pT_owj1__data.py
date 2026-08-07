
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.100000e+02, 4.800000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.100000e+02, 4.800000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.000000e+01, 7.000000e+01, 1.025000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.100000e+02, 4.800000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   2.000000e+02, 2.600000e+02, 3.600000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   2.000000e+02, 2.600000e+02, 3.600000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 5.500000e+01, 8.500000e+01, 1.200000e+02, 1.600000e+02,
                                                                                   2.000000e+02, 2.600000e+02, 3.600000e+02, 6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.868141e-01, 2.707826e-02, 4.871546e-03, 1.560080e-03, 8.841705e-04,
                                                                                   5.627503e-04, 3.077597e-04, 1.206024e-04],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.848965e-01, 2.664804e-02, 4.716799e-03, 1.474478e-03, 8.277796e-04,
                                                                                   5.365362e-04, 2.970622e-04, 1.168324e-04],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.788701e-01, 2.616187e-02, 4.637202e-03, 1.466766e-03, 8.187390e-04,
                                                                                   5.259233e-04, 2.857108e-04, 1.139391e-04],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                     [1.500000e+01, 1.500000e+01, 1.750000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      3.000000e+01, 5.000000e+01, 1.200000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.957240e-05, 2.273638e-05, 8.955936e-06, 4.741697e-06, 3.562947e-06,
                                                                                      2.319127e-06, 1.325634e-06, 5.347838e-07],
                                                                                     [5.957240e-05, 2.273638e-05, 8.955936e-06, 4.741697e-06, 3.562947e-06,
                                                                                      2.319127e-06, 1.325634e-06, 5.347838e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.157480e-04, 4.406545e-05, 1.721217e-05, 9.003204e-06, 6.740181e-06,
                                                                                      4.414440e-06, 2.544923e-06, 1.027783e-06],
                                                                                     [1.157480e-04, 4.406545e-05, 1.721217e-05, 9.003204e-06, 6.740181e-06,
                                                                                      4.414440e-06, 2.544923e-06, 1.027783e-06],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.092537e-04, 4.191296e-05, 1.639867e-05, 8.621030e-06, 6.424553e-06,
                                                                                      4.208126e-06, 2.394022e-06, 9.747227e-07],
                                                                                     [1.092537e-04, 4.191296e-05, 1.639867e-05, 8.621030e-06, 6.424553e-06,
                                                                                      4.208126e-06, 2.394022e-06, 9.747227e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.897353e-01, 9.841120e-01, 9.682345e-01, 9.451297e-01, 9.362217e-01,
                                                                                   9.534179e-01, 9.652407e-01, 9.687403e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.574764e-01, 9.661577e-01, 9.518954e-01, 9.401864e-01, 9.259967e-01,
                                                                                   9.345589e-01, 9.283568e-01, 9.447499e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.188860e-04, 8.396544e-04, 1.838418e-03, 3.039393e-03, 4.029706e-03,
                                                                                      4.121059e-03, 4.307367e-03, 4.434272e-03],
                                                                                     [3.188860e-04, 8.396544e-04, 1.838418e-03, 3.039393e-03, 4.029706e-03,
                                                                                      4.121059e-03, 4.307367e-03, 4.434272e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.195892e-04, 1.627337e-03, 3.533205e-03, 5.770989e-03, 7.623169e-03,
                                                                                      7.844403e-03, 8.269189e-03, 8.522078e-03],
                                                                                     [6.195892e-04, 1.627337e-03, 3.533205e-03, 5.770989e-03, 7.623169e-03,
                                                                                      7.844403e-03, 8.269189e-03, 8.522078e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.848258e-04, 1.547845e-03, 3.366215e-03, 5.526018e-03, 7.266192e-03,
                                                                                      7.477785e-03, 7.778868e-03, 8.082117e-03],
                                                                                     [5.848258e-04, 1.547845e-03, 3.366215e-03, 5.526018e-03, 7.266192e-03,
                                                                                      7.477785e-03, 7.778868e-03, 8.082117e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
