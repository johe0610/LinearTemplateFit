
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.600000e+02,
                                                                                   3.000000e+02, 3.450000e+02, 4.000000e+02, 4.750000e+02, 6.000000e+02,
                                                                                   7.900000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.600000e+02,
                                                                                   3.000000e+02, 3.450000e+02, 4.000000e+02, 4.750000e+02, 6.000000e+02,
                                                                                   7.900000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.600000e+02,
                                                                                   3.000000e+02, 3.450000e+02, 4.000000e+02, 4.750000e+02, 6.000000e+02,
                                                                                   7.900000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   2.800000e+02, 3.200000e+02, 3.700000e+02, 4.300000e+02, 5.200000e+02,
                                                                                   6.800000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   2.800000e+02, 3.200000e+02, 3.700000e+02, 4.300000e+02, 5.200000e+02,
                                                                                   6.800000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   2.800000e+02, 3.200000e+02, 3.700000e+02, 4.300000e+02, 5.200000e+02,
                                                                                   6.800000e+02, 9.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.751919e-01, 1.095298e-01, 5.808391e-02, 2.977021e-02, 1.545289e-02,
                                                                                   8.373120e-03, 4.648459e-03, 2.649389e-03, 1.378131e-03, 4.987653e-04,
                                                                                   1.210696e-04],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.748101e-01, 1.090039e-01, 5.731902e-02, 2.923055e-02, 1.508623e-02,
                                                                                   8.151301e-03, 4.471867e-03, 2.534355e-03, 1.317276e-03, 4.654682e-04,
                                                                                   1.108808e-04],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.658942e-01, 1.054863e-01, 5.653599e-02, 2.902245e-02, 1.497760e-02,
                                                                                   8.082206e-03, 4.426772e-03, 2.497669e-03, 1.290646e-03, 4.551130e-04,
                                                                                   1.071334e-04],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01,
                                                                                      2.000000e+01, 2.500000e+01, 3.000000e+01, 4.500000e+01, 8.000000e+01,
                                                                                      1.100000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.984826e-05, 3.946601e-05, 2.879569e-05, 2.066125e-05, 1.492451e-05,
                                                                                      1.101353e-05, 7.357579e-06, 5.083497e-06, 2.999994e-06, 1.356099e-06,
                                                                                      5.681102e-07],
                                                                                     [4.984826e-05, 3.946601e-05, 2.879569e-05, 2.066125e-05, 1.492451e-05,
                                                                                      1.101353e-05, 7.357579e-06, 5.083497e-06, 2.999994e-06, 1.356099e-06,
                                                                                      5.681102e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.724117e-05, 7.690019e-05, 5.586912e-05, 3.999120e-05, 2.880330e-05,
                                                                                      2.122977e-05, 1.411266e-05, 9.712793e-06, 5.739922e-06, 2.562964e-06,
                                                                                      1.062126e-06],
                                                                                     [9.724117e-05, 7.690019e-05, 5.586912e-05, 3.999120e-05, 2.880330e-05,
                                                                                      2.122977e-05, 1.411266e-05, 9.712793e-06, 5.739922e-06, 2.562964e-06,
                                                                                      1.062126e-06],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.091285e-05, 7.259528e-05, 5.324348e-05, 3.823501e-05, 2.753747e-05,
                                                                                      2.028730e-05, 1.347971e-05, 9.252187e-06, 5.451338e-06, 2.433705e-06,
                                                                                      1.004945e-06],
                                                                                     [9.091285e-05, 7.259528e-05, 5.324348e-05, 3.823501e-05, 2.753747e-05,
                                                                                      2.028730e-05, 1.347971e-05, 9.252187e-06, 5.451338e-06, 2.433705e-06,
                                                                                      1.004945e-06],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.978207e-01, 9.951986e-01, 9.868313e-01, 9.818725e-01, 9.762724e-01,
                                                                                   9.735082e-01, 9.620106e-01, 9.565809e-01, 9.558424e-01, 9.332409e-01,
                                                                                   9.158434e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.469285e-01, 9.630831e-01, 9.733503e-01, 9.748823e-01, 9.692426e-01,
                                                                                   9.652562e-01, 9.523096e-01, 9.427340e-01, 9.365191e-01, 9.124793e-01,
                                                                                   8.848910e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.845352e-04, 3.603221e-04, 4.957602e-04, 6.940243e-04, 9.658070e-04,
                                                                                      1.315344e-03, 1.582800e-03, 1.918743e-03, 2.176857e-03, 2.718912e-03,
                                                                                      4.692427e-03],
                                                                                     [2.845352e-04, 3.603221e-04, 4.957602e-04, 6.940243e-04, 9.658070e-04,
                                                                                      1.315344e-03, 1.582800e-03, 1.918743e-03, 2.176857e-03, 2.718912e-03,
                                                                                      4.692427e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.550552e-04, 7.020938e-04, 9.618691e-04, 1.343329e-03, 1.863943e-03,
                                                                                      2.535467e-03, 3.035987e-03, 3.666050e-03, 4.165005e-03, 5.138617e-03,
                                                                                      8.772855e-03],
                                                                                     [5.550552e-04, 7.020938e-04, 9.618691e-04, 1.343329e-03, 1.863943e-03,
                                                                                      2.535467e-03, 3.035987e-03, 3.666050e-03, 4.165005e-03, 5.138617e-03,
                                                                                      8.772855e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.189330e-04, 6.627902e-04, 9.166649e-04, 1.284338e-03, 1.782027e-03,
                                                                                      2.422908e-03, 2.899823e-03, 3.492197e-03, 3.955602e-03, 4.879459e-03,
                                                                                      8.300556e-03],
                                                                                     [5.189330e-04, 6.627902e-04, 9.166649e-04, 1.284338e-03, 1.782027e-03,
                                                                                      2.422908e-03, 2.899823e-03, 3.492197e-03, 3.955602e-03, 4.879459e-03,
                                                                                      8.300556e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
