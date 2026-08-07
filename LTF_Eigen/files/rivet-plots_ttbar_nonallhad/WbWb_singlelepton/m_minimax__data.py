
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.700000e+02,
                                                                                   3.300000e+02, 4.000000e+02, 4.800000e+02, 5.700000e+02, 6.700000e+02,
                                                                                   8.100000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.700000e+02,
                                                                                   3.300000e+02, 4.000000e+02, 4.800000e+02, 5.700000e+02, 6.700000e+02,
                                                                                   8.100000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.700000e+02,
                                                                                   3.300000e+02, 4.000000e+02, 4.800000e+02, 5.700000e+02, 6.700000e+02,
                                                                                   8.100000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02, 6.200000e+02,
                                                                                   7.200000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02, 6.200000e+02,
                                                                                   7.200000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02, 6.200000e+02,
                                                                                   7.200000e+02, 9.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.345153e-02, 1.014738e-01, 2.508348e-01, 2.233188e-02, 7.762231e-03,
                                                                                   3.212975e-03, 1.600102e-03, 8.477193e-04, 4.645603e-04, 2.512782e-04,
                                                                                   1.170337e-04],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.347888e-02, 1.087931e-01, 2.449570e-01, 2.042070e-02, 6.849521e-03,
                                                                                   2.751221e-03, 1.364203e-03, 7.310053e-04, 3.992357e-04, 2.213370e-04,
                                                                                   1.033316e-04],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.230944e-02, 9.059943e-02, 2.502052e-01, 2.097202e-02, 6.857762e-03,
                                                                                   2.705209e-03, 1.324601e-03, 7.049995e-04, 3.905734e-04, 2.073947e-04,
                                                                                   9.968594e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.381383e-05, 3.796843e-05, 5.983674e-05, 1.782153e-05, 8.579537e-06,
                                                                                      5.520441e-06, 3.375587e-06, 2.457515e-06, 1.627988e-06, 1.198278e-06,
                                                                                      6.093593e-07],
                                                                                     [1.381383e-05, 3.796843e-05, 5.983674e-05, 1.782153e-05, 8.579537e-06,
                                                                                      5.520441e-06, 3.375587e-06, 2.457515e-06, 1.627988e-06, 1.198278e-06,
                                                                                      6.093593e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.700489e-05, 7.678733e-05, 1.154948e-04, 3.327814e-05, 1.573870e-05,
                                                                                      9.976421e-06, 6.088616e-06, 4.453149e-06, 2.947042e-06, 2.197294e-06,
                                                                                      1.118660e-06],
                                                                                     [2.700489e-05, 7.678733e-05, 1.154948e-04, 3.327814e-05, 1.573870e-05,
                                                                                      9.976421e-06, 6.088616e-06, 4.453149e-06, 2.947042e-06, 2.197294e-06,
                                                                                      1.118660e-06],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.476504e-05, 6.724379e-05, 1.120081e-04, 3.237369e-05, 1.511468e-05,
                                                                                      9.496755e-06, 5.754211e-06, 4.198502e-06, 2.797172e-06, 2.040591e-06,
                                                                                      1.056296e-06],
                                                                                     [2.476504e-05, 6.724379e-05, 1.120081e-04, 3.237369e-05, 1.511468e-05,
                                                                                      9.496755e-06, 5.754211e-06, 4.198502e-06, 2.797172e-06, 2.040591e-06,
                                                                                      1.056296e-06],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.002033e+00, 1.072130e+00, 9.765670e-01, 9.144192e-01, 8.824165e-01,
                                                                                   8.562846e-01, 8.525725e-01, 8.623200e-01, 8.593840e-01, 8.808444e-01,
                                                                                   8.829218e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.150959e-01, 8.928357e-01, 9.974900e-01, 9.391068e-01, 8.834782e-01,
                                                                                   8.419639e-01, 8.278229e-01, 8.316426e-01, 8.407378e-01, 8.253589e-01,
                                                                                   8.517712e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.026934e-03, 3.741698e-04, 2.385504e-04, 7.980309e-04, 1.105293e-03,
                                                                                      1.718171e-03, 2.109607e-03, 2.898973e-03, 3.504363e-03, 4.768730e-03,
                                                                                      5.206699e-03],
                                                                                     [1.026934e-03, 3.741698e-04, 2.385504e-04, 7.980309e-04, 1.105293e-03,
                                                                                      1.718171e-03, 2.109607e-03, 2.898973e-03, 3.504363e-03, 4.768730e-03,
                                                                                      5.206699e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.007570e-03, 7.567207e-04, 4.604417e-04, 1.490163e-03, 2.027600e-03,
                                                                                      3.105042e-03, 3.805142e-03, 5.253094e-03, 6.343723e-03, 8.744467e-03,
                                                                                      9.558443e-03],
                                                                                     [2.007570e-03, 7.567207e-04, 4.604417e-04, 1.490163e-03, 2.027600e-03,
                                                                                      3.105042e-03, 3.805142e-03, 5.253094e-03, 6.343723e-03, 8.744467e-03,
                                                                                      9.558443e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.841057e-03, 6.626714e-04, 4.465413e-04, 1.449663e-03, 1.947208e-03,
                                                                                      2.955751e-03, 3.596153e-03, 4.952703e-03, 6.021117e-03, 8.120844e-03,
                                                                                      9.025571e-03],
                                                                                     [1.841057e-03, 6.626714e-04, 4.465413e-04, 1.449663e-03, 1.947208e-03,
                                                                                      2.955751e-03, 3.596153e-03, 4.952703e-03, 6.021117e-03, 8.120844e-03,
                                                                                      9.025571e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
