
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.700000e+02, 2.700000e+02,
                                                                                   4.700000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.700000e+02, 2.700000e+02,
                                                                                   4.700000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.700000e+02, 2.700000e+02,
                                                                                   4.700000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 2.000000e+02,
                                                                                   3.400000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 2.000000e+02,
                                                                                   3.400000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 2.000000e+02,
                                                                                   3.400000e+02, 6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.102170e-01, 9.759701e-02, 2.801258e-02, 7.624855e-03, 9.801212e-04,
                                                                                   3.709705e-05],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.106333e-01, 9.514877e-02, 2.713737e-02, 7.331318e-03, 9.429622e-04,
                                                                                   3.501899e-05],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.958413e-01, 9.414795e-02, 2.705482e-02, 7.406912e-03, 9.620479e-04,
                                                                                   3.535256e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01, 7.000000e+01,
                                                                                      1.300000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.098129e-05, 3.731655e-05, 2.006999e-05, 8.588084e-06, 2.031329e-06,
                                                                                      2.960905e-07],
                                                                                     [7.098129e-05, 3.731655e-05, 2.006999e-05, 8.588084e-06, 2.031329e-06,
                                                                                      2.960905e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.387204e-04, 7.197416e-05, 3.858092e-05, 1.645585e-05, 3.890688e-06,
                                                                                      5.644109e-07],
                                                                                     [1.387204e-04, 7.197416e-05, 3.858092e-05, 1.645585e-05, 3.890688e-06,
                                                                                      5.644109e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.299286e-04, 6.869549e-05, 3.696396e-05, 1.586161e-05, 3.772746e-06,
                                                                                      5.440039e-07],
                                                                                     [1.299286e-04, 6.869549e-05, 3.696396e-05, 1.586161e-05, 3.772746e-06,
                                                                                      5.440039e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.001342e+00, 9.749148e-01, 9.687565e-01, 9.615026e-01, 9.620873e-01,
                                                                                   9.439831e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.536592e-01, 9.646602e-01, 9.658096e-01, 9.714168e-01, 9.815601e-01,
                                                                                   9.529750e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.288117e-04, 3.823534e-04, 7.164635e-04, 1.126328e-03, 2.072528e-03,
                                                                                      7.981511e-03],
                                                                                     [2.288117e-04, 3.823534e-04, 7.164635e-04, 1.126328e-03, 2.072528e-03,
                                                                                      7.981511e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.471721e-04, 7.374628e-04, 1.377271e-03, 2.158185e-03, 3.969599e-03,
                                                                                      1.521444e-02],
                                                                                     [4.471721e-04, 7.374628e-04, 1.377271e-03, 2.158185e-03, 3.969599e-03,
                                                                                      1.521444e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.188313e-04, 7.038688e-04, 1.319549e-03, 2.080251e-03, 3.849265e-03,
                                                                                      1.466434e-02],
                                                                                     [4.188313e-04, 7.038688e-04, 1.319549e-03, 2.080251e-03, 3.849265e-03,
                                                                                      1.466434e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
