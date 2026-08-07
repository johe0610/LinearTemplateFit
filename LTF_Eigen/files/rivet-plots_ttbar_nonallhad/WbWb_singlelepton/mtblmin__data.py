
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 2.000000e+02, 2.800000e+02,
                                                                                   3.800000e+02, 5.050000e+02, 7.350000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 2.000000e+02, 2.800000e+02,
                                                                                   3.800000e+02, 5.050000e+02, 7.350000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 7.500000e+01, 1.300000e+02, 2.000000e+02, 2.800000e+02,
                                                                                   3.800000e+02, 5.050000e+02, 7.350000e+02, 9.500000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.400000e+02,
                                                                                   3.200000e+02, 4.400000e+02, 5.700000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.400000e+02,
                                                                                   3.200000e+02, 4.400000e+02, 5.700000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 5.000000e+01, 1.000000e+02, 1.600000e+02, 2.400000e+02,
                                                                                   3.200000e+02, 4.400000e+02, 5.700000e+02, 9.000000e+02, 1.000000e+03],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [7.259052e-03, 8.265880e-02, 1.738610e-01, 1.838656e-02, 8.506876e-04,
                                                                                   1.376883e-04, 2.199800e-05, 2.291670e-06, 2.821682e-07],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [7.236387e-03, 8.416467e-02, 1.746135e-01, 1.564499e-02, 4.813168e-04,
                                                                                   5.850990e-05, 6.903850e-06, 3.914045e-07, 2.153608e-08],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.955971e-03, 7.905472e-02, 1.669213e-01, 1.784565e-02, 4.903266e-04,
                                                                                   6.184058e-05, 7.154208e-06, 5.228449e-07, 3.965292e-08],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                     [2.500000e+01, 2.500000e+01, 3.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      6.000000e+01, 6.500000e+01, 1.650000e+02, 5.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.105031e-06, 3.068008e-05, 4.063466e-05, 1.146596e-05, 2.474672e-06,
                                                                                      8.178820e-07, 3.175779e-07, 6.534689e-08, 4.069618e-08],
                                                                                     [9.105031e-06, 3.068008e-05, 4.063466e-05, 1.146596e-05, 2.474672e-06,
                                                                                      8.178820e-07, 3.175779e-07, 6.534689e-08, 4.069618e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.776151e-05, 6.046705e-05, 7.953390e-05, 2.066080e-05, 3.654531e-06,
                                                                                      1.046811e-06, 3.476742e-07, 5.299577e-08, 2.153608e-08],
                                                                                     [1.776151e-05, 6.046705e-05, 7.953390e-05, 2.066080e-05, 3.654531e-06,
                                                                                      1.046811e-06, 3.476742e-07, 5.299577e-08, 2.153608e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.671072e-05, 5.623921e-05, 7.462868e-05, 2.117364e-05, 3.529635e-06,
                                                                                      1.034699e-06, 3.421083e-07, 5.795638e-08, 2.803885e-08],
                                                                                     [1.671072e-05, 5.623921e-05, 7.462868e-05, 2.117364e-05, 3.529635e-06,
                                                                                      1.034699e-06, 3.421083e-07, 5.795638e-08, 2.803885e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.968777e-01, 1.018218e+00, 1.004328e+00, 8.508927e-01, 5.657974e-01,
                                                                                   4.249446e-01, 3.138399e-01, 1.707944e-01, 7.632355e-02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.582479e-01, 9.563981e-01, 9.600848e-01, 9.705812e-01, 5.763886e-01,
                                                                                   4.491346e-01, 3.252208e-01, 2.281502e-01, 1.405294e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.254300e-03, 3.711653e-04, 2.337192e-04, 6.236055e-04, 2.909026e-03,
                                                                                      5.940098e-03, 1.443667e-02, 2.851497e-02, 1.442267e-01],
                                                                                     [1.254300e-03, 3.711653e-04, 2.337192e-04, 6.236055e-04, 2.909026e-03,
                                                                                      5.940098e-03, 1.443667e-02, 2.851497e-02, 1.442267e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.446808e-03, 7.315259e-04, 4.574568e-04, 1.123690e-03, 4.295973e-03,
                                                                                      7.602759e-03, 1.580481e-02, 2.312539e-02, 7.632355e-02],
                                                                                     [2.446808e-03, 7.315259e-04, 4.574568e-04, 1.123690e-03, 4.295973e-03,
                                                                                      7.602759e-03, 1.580481e-02, 2.312539e-02, 7.632355e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.302053e-03, 6.803778e-04, 4.292434e-04, 1.151582e-03, 4.149155e-03,
                                                                                      7.514792e-03, 1.555179e-02, 2.529002e-02, 9.936928e-02],
                                                                                     [2.302053e-03, 6.803778e-04, 4.292434e-04, 1.151582e-03, 4.149155e-03,
                                                                                      7.514792e-03, 1.555179e-02, 2.529002e-02, 9.936928e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
