
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 5.000000e+01, 7.000000e+01, 9.000000e+01, 1.100000e+02,
                                                                                   1.350000e+02, 1.750000e+02, 2.750000e+02, 4.750000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 5.000000e+01, 7.000000e+01, 9.000000e+01, 1.100000e+02,
                                                                                   1.350000e+02, 1.750000e+02, 2.750000e+02, 4.750000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 5.000000e+01, 7.000000e+01, 9.000000e+01, 1.100000e+02,
                                                                                   1.350000e+02, 1.750000e+02, 2.750000e+02, 4.750000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 6.000000e+01, 8.000000e+01, 1.000000e+02,
                                                                                   1.200000e+02, 1.500000e+02, 2.000000e+02, 3.500000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 6.000000e+01, 8.000000e+01, 1.000000e+02,
                                                                                   1.200000e+02, 1.500000e+02, 2.000000e+02, 3.500000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 6.000000e+01, 8.000000e+01, 1.000000e+02,
                                                                                   1.200000e+02, 1.500000e+02, 2.000000e+02, 3.500000e+02, 6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [7.917461e-04, 2.593568e-03, 3.875329e-03, 4.456527e-03, 3.999366e-03,
                                                                                   2.078903e-03, 8.506415e-05, 4.466084e-06, 2.061225e-07],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.334528e-04, 2.721325e-03, 4.040781e-03, 4.625922e-03, 4.070429e-03,
                                                                                   2.010041e-03, 7.488555e-05, 4.440175e-06, 2.196323e-07],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [7.536046e-04, 2.483866e-03, 3.738835e-03, 4.353556e-03, 3.970042e-03,
                                                                                   2.161318e-03, 9.977425e-05, 4.427578e-06, 2.023371e-07],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.670766e-06, 4.277383e-06, 5.229846e-06, 5.609107e-06, 5.315088e-06,
                                                                                      3.130805e-06, 4.907700e-07, 6.487973e-08, 1.086111e-08],
                                                                                     [1.670766e-06, 4.277383e-06, 5.229846e-06, 5.609107e-06, 5.315088e-06,
                                                                                      3.130805e-06, 4.907700e-07, 6.487973e-08, 1.086111e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.662449e-06, 1.191663e-05, 1.452374e-05, 1.553902e-05, 1.458307e-05,
                                                                                      8.376766e-06, 1.251004e-06, 1.759182e-07, 3.016882e-08],
                                                                                     [4.662449e-06, 1.191663e-05, 1.452374e-05, 1.553902e-05, 1.458307e-05,
                                                                                      8.376766e-06, 1.251004e-06, 1.759182e-07, 3.016882e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.252629e-06, 1.092586e-05, 1.340615e-05, 1.446627e-05, 1.382166e-05,
                                                                                      8.329614e-06, 1.390162e-06, 1.685492e-07, 2.779315e-08],
                                                                                     [4.252629e-06, 1.092586e-05, 1.340615e-05, 1.446627e-05, 1.382166e-05,
                                                                                      8.329614e-06, 1.390162e-06, 1.685492e-07, 2.779315e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.052677e+00, 1.049259e+00, 1.042694e+00, 1.038011e+00, 1.017769e+00,
                                                                                   9.668758e-01, 8.803421e-01, 9.941987e-01, 1.065543e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.518261e-01, 9.577023e-01, 9.647787e-01, 9.768943e-01, 9.926678e-01,
                                                                                   1.039644e+00, 1.172929e+00, 9.913781e-01, 9.816352e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.110230e-03, 1.649227e-03, 1.349523e-03, 1.258627e-03, 1.328983e-03,
                                                                                      1.505989e-03, 5.769410e-03, 1.452721e-02, 5.269250e-02],
                                                                                     [2.110230e-03, 1.649227e-03, 1.349523e-03, 1.258627e-03, 1.328983e-03,
                                                                                      1.505989e-03, 5.769410e-03, 1.452721e-02, 5.269250e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.888818e-03, 4.594686e-03, 3.747744e-03, 3.486800e-03, 3.646345e-03,
                                                                                      4.029416e-03, 1.470659e-02, 3.938981e-02, 1.463635e-01],
                                                                                     [5.888818e-03, 4.594686e-03, 3.747744e-03, 3.486800e-03, 3.646345e-03,
                                                                                      4.029416e-03, 1.470659e-02, 3.938981e-02, 1.463635e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.371203e-03, 4.212675e-03, 3.459358e-03, 3.246086e-03, 3.455963e-03,
                                                                                      4.006735e-03, 1.634251e-02, 3.773982e-02, 1.348380e-01],
                                                                                     [5.371203e-03, 4.212675e-03, 3.459358e-03, 3.246086e-03, 3.455963e-03,
                                                                                      4.006735e-03, 1.634251e-02, 3.773982e-02, 1.348380e-01],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
