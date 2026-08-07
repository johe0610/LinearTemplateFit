
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
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.658043e+00, 1.126435e+01, 1.344936e+01, 9.466085e+00, 2.935235e+00,
                                                                                   4.320178e-01, 4.382769e-02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.636008e+00, 1.132475e+01, 1.344979e+01, 9.385106e+00, 2.853415e+00,
                                                                                   4.072163e-01, 4.117741e-02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.525740e+00, 1.069865e+01, 1.285957e+01, 9.121779e+00, 2.831838e+00,
                                                                                   4.123935e-01, 4.119958e-02],
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
                                                                                     [1.950361e-03, 4.624992e-03, 5.050864e-03, 3.282713e-03, 1.364206e-03,
                                                                                      4.547729e-04, 1.031100e-04],
                                                                                     [1.950361e-03, 4.624992e-03, 5.050864e-03, 3.282713e-03, 1.364206e-03,
                                                                                      4.547729e-04, 1.031100e-04],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.794049e-03, 9.056658e-03, 9.865113e-03, 6.384406e-03, 2.626813e-03,
                                                                                      8.622529e-04, 1.952534e-04],
                                                                                     [3.794049e-03, 9.056658e-03, 9.865113e-03, 6.384406e-03, 2.626813e-03,
                                                                                      8.622529e-04, 1.952534e-04],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.564971e-03, 8.448818e-03, 9.257552e-03, 6.040401e-03, 2.511060e-03,
                                                                                      8.324753e-04, 1.873946e-04],
                                                                                     [3.564971e-03, 8.448818e-03, 9.257552e-03, 6.040401e-03, 2.511060e-03,
                                                                                      8.324753e-04, 1.873946e-04],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.917101e-01, 1.005362e+00, 1.000032e+00, 9.914454e-01, 9.721249e-01,
                                                                                   9.425915e-01, 9.395296e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.502254e-01, 9.497796e-01, 9.561474e-01, 9.636274e-01, 9.647739e-01,
                                                                                   9.545753e-01, 9.400354e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.337583e-04, 4.105867e-04, 3.755468e-04, 3.467868e-04, 4.647689e-04,
                                                                                      1.052672e-03, 2.352622e-03],
                                                                                     [7.337583e-04, 4.105867e-04, 3.755468e-04, 3.467868e-04, 4.647689e-04,
                                                                                      1.052672e-03, 2.352622e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.427384e-03, 8.040107e-04, 7.335006e-04, 6.744505e-04, 8.949243e-04,
                                                                                      1.995874e-03, 4.455024e-03],
                                                                                     [1.427384e-03, 8.040107e-04, 7.335006e-04, 6.744505e-04, 8.949243e-04,
                                                                                      1.995874e-03, 4.455024e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.341201e-03, 7.500493e-04, 6.883266e-04, 6.381097e-04, 8.554886e-04,
                                                                                      1.926947e-03, 4.275712e-03],
                                                                                     [1.341201e-03, 7.500493e-04, 6.883266e-04, 6.381097e-04, 8.554886e-04,
                                                                                      1.926947e-03, 4.275712e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
