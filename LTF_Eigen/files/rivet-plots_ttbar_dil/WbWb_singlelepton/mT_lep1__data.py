
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 6.500000e+01, 1.150000e+02, 1.700000e+02, 2.500000e+02,
                                                                                   4.000000e+02, 7.000000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 6.500000e+01, 1.150000e+02, 1.700000e+02, 2.500000e+02,
                                                                                   4.000000e+02, 7.000000e+02, 9.500000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 6.500000e+01, 1.150000e+02, 1.700000e+02, 2.500000e+02,
                                                                                   4.000000e+02, 7.000000e+02, 9.500000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 9.000000e+01, 1.400000e+02, 2.000000e+02,
                                                                                   3.000000e+02, 5.000000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 9.000000e+01, 1.400000e+02, 2.000000e+02,
                                                                                   3.000000e+02, 5.000000e+02, 9.000000e+02, 1.000000e+03],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 9.000000e+01, 1.400000e+02, 2.000000e+02,
                                                                                   3.000000e+02, 5.000000e+02, 9.000000e+02, 1.000000e+03],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.516172e-03, 2.868028e-03, 2.266643e-03, 8.792827e-04, 2.161848e-04,
                                                                                   2.624470e-05, 1.216181e-06, 9.102146e-08],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.546325e-03, 2.949045e-03, 2.330370e-03, 8.950525e-04, 2.185340e-04,
                                                                                   2.680483e-05, 1.286964e-06, 7.248633e-08],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.507501e-03, 2.815171e-03, 2.223529e-03, 8.663961e-04, 2.164701e-04,
                                                                                   2.600272e-05, 1.333501e-06, 9.542974e-08],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                     [2.000000e+01, 2.500000e+01, 2.500000e+01, 3.000000e+01, 5.000000e+01,
                                                                                      1.000000e+02, 2.000000e+02, 5.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.312459e-06, 2.843681e-06, 2.528705e-06, 1.440045e-06, 5.546276e-07,
                                                                                      1.378228e-07, 2.151777e-08, 1.196313e-08],
                                                                                     [2.312459e-06, 2.843681e-06, 2.528705e-06, 1.440045e-06, 5.546276e-07,
                                                                                      1.378228e-07, 2.151777e-08, 1.196313e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.352641e-06, 7.840589e-06, 6.973131e-06, 3.950740e-06, 1.517975e-06,
                                                                                      3.797077e-07, 6.078388e-08, 3.106347e-08],
                                                                                     [6.352641e-06, 7.840589e-06, 6.973131e-06, 3.950740e-06, 1.517975e-06,
                                                                                      3.797077e-07, 6.078388e-08, 3.106347e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.014900e-06, 7.352434e-06, 6.537343e-06, 3.731115e-06, 1.448407e-06,
                                                                                      3.579949e-07, 5.887074e-08, 3.017754e-08],
                                                                                     [6.014900e-06, 7.352434e-06, 6.537343e-06, 3.731115e-06, 1.448407e-06,
                                                                                      3.579949e-07, 5.887074e-08, 3.017754e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.019888e+00, 1.028248e+00, 1.028115e+00, 1.017935e+00, 1.010867e+00,
                                                                                   1.021343e+00, 1.058201e+00, 7.963653e-01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.942810e-01, 9.815703e-01, 9.809789e-01, 9.853442e-01, 1.001320e+00,
                                                                                   9.907799e-01, 1.096466e+00, 1.048431e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.525196e-03, 9.915109e-04, 1.115617e-03, 1.637750e-03, 2.565525e-03,
                                                                                      5.251453e-03, 1.769290e-02, 1.314320e-01],
                                                                                     [1.525196e-03, 9.915109e-04, 1.115617e-03, 1.637750e-03, 2.565525e-03,
                                                                                      5.251453e-03, 1.769290e-02, 1.314320e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.189921e-03, 2.733791e-03, 3.076413e-03, 4.493140e-03, 7.021655e-03,
                                                                                      1.446798e-02, 4.997930e-02, 3.412763e-01],
                                                                                     [4.189921e-03, 2.733791e-03, 3.076413e-03, 4.493140e-03, 7.021655e-03,
                                                                                      1.446798e-02, 4.997930e-02, 3.412763e-01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.967162e-03, 2.563585e-03, 2.884152e-03, 4.243362e-03, 6.699856e-03,
                                                                                      1.364066e-02, 4.840623e-02, 3.315431e-01],
                                                                                     [3.967162e-03, 2.563585e-03, 2.884152e-03, 4.243362e-03, 6.699856e-03,
                                                                                      1.364066e-02, 4.840623e-02, 3.315431e-01],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
