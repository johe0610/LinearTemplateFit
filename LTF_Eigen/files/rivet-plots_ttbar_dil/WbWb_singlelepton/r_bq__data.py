
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e-01, 5.500000e-01, 8.500000e-01, 1.250000e+00, 1.950000e+00,
                                                                                   3.000000e+00, 4.800000e+00],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e-01, 7.000000e-01, 1.000000e+00, 1.500000e+00,
                                                                                   2.400000e+00, 3.600000e+00, 6.000000e+00],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.304139e-02, 1.796740e-01, 2.617286e-01, 2.330064e-01, 1.079529e-01,
                                                                                   2.468659e-02, 3.244371e-03],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.453245e-02, 1.873077e-01, 2.723402e-01, 2.391349e-01, 1.089612e-01,
                                                                                   2.449675e-02, 3.221625e-03],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.189989e-02, 1.737342e-01, 2.553443e-01, 2.285796e-01, 1.080190e-01,
                                                                                   2.489769e-02, 3.288526e-03],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                     [2.000000e-01, 1.500000e-01, 1.500000e-01, 2.500000e-01, 4.500000e-01,
                                                                                      6.000000e-01, 1.200000e+00],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.078024e-04, 2.903085e-04, 3.505250e-04, 2.563220e-04, 1.302327e-04,
                                                                                      5.414872e-05, 1.398414e-05],
                                                                                     [1.078024e-04, 2.903085e-04, 3.505250e-04, 2.563220e-04, 1.302327e-04,
                                                                                      5.414872e-05, 1.398414e-05],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.993660e-04, 8.060212e-04, 9.725459e-04, 7.062449e-04, 3.558582e-04,
                                                                                      1.467058e-04, 3.801461e-05],
                                                                                     [2.993660e-04, 8.060212e-04, 9.725459e-04, 7.062449e-04, 3.558582e-04,
                                                                                      1.467058e-04, 3.801461e-05],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.763147e-04, 7.450497e-04, 9.036863e-04, 6.626414e-04, 3.398452e-04,
                                                                                      1.418400e-04, 3.682233e-05],
                                                                                     [2.763147e-04, 7.450497e-04, 9.036863e-04, 6.626414e-04, 3.398452e-04,
                                                                                      1.418400e-04, 3.682233e-05],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.045127e+00, 1.042486e+00, 1.040544e+00, 1.026302e+00, 1.009340e+00,
                                                                                   9.923100e-01, 9.929891e-01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.654524e-01, 9.669412e-01, 9.756072e-01, 9.810014e-01, 1.000612e+00,
                                                                                   1.008551e+00, 1.013610e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.262647e-03, 1.615751e-03, 1.339269e-03, 1.100064e-03, 1.206384e-03,
                                                                                      2.193447e-03, 4.310278e-03],
                                                                                     [3.262647e-03, 1.615751e-03, 1.339269e-03, 1.100064e-03, 1.206384e-03,
                                                                                      2.193447e-03, 4.310278e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.060333e-03, 4.486020e-03, 3.715856e-03, 3.031011e-03, 3.296421e-03,
                                                                                      5.942732e-03, 1.171710e-02],
                                                                                     [9.060333e-03, 4.486020e-03, 3.715856e-03, 3.031011e-03, 3.296421e-03,
                                                                                      5.942732e-03, 1.171710e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.362684e-03, 4.146675e-03, 3.452761e-03, 2.843876e-03, 3.148088e-03,
                                                                                      5.745630e-03, 1.134961e-02],
                                                                                     [8.362684e-03, 4.146675e-03, 3.452761e-03, 2.843876e-03, 3.148088e-03,
                                                                                      5.745630e-03, 1.134961e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
