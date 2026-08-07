
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 5.000000e+01, 7.000000e+01, 9.000000e+01, 1.100000e+02,
                                                                                   1.350000e+02, 1.750000e+02, 2.750000e+02, 4.750000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 5.000000e+01, 7.000000e+01, 9.000000e+01, 1.100000e+02,
                                                                                   1.350000e+02, 1.750000e+02, 2.750000e+02, 4.750000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+01, 5.000000e+01, 7.000000e+01, 9.000000e+01, 1.100000e+02,
                                                                                   1.350000e+02, 1.750000e+02, 2.750000e+02, 4.750000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 6.000000e+01, 8.000000e+01, 1.000000e+02,
                                                                                   1.200000e+02, 1.500000e+02, 2.000000e+02, 3.500000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 6.000000e+01, 8.000000e+01, 1.000000e+02,
                                                                                   1.200000e+02, 1.500000e+02, 2.000000e+02, 3.500000e+02, 6.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 4.000000e+01, 6.000000e+01, 8.000000e+01, 1.000000e+02,
                                                                                   1.200000e+02, 1.500000e+02, 2.000000e+02, 3.500000e+02, 6.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.522975e-02, 1.135648e-01, 1.687738e-01, 1.881216e-01, 1.592099e-01,
                                                                                   7.732194e-02, 2.865862e-03, 1.317258e-04, 5.363849e-06],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.585466e-02, 1.152131e-01, 1.706269e-01, 1.882243e-01, 1.567284e-01,
                                                                                   7.190591e-02, 2.353667e-03, 1.253431e-04, 4.993784e-06],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.295560e-02, 1.062712e-01, 1.591887e-01, 1.787061e-01, 1.530803e-01,
                                                                                   7.801931e-02, 3.329485e-03, 1.301973e-04, 5.497213e-06],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                     [2.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01, 1.000000e+01,
                                                                                      1.500000e+01, 2.500000e+01, 7.500000e+01, 1.250000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.237953e-05, 5.685410e-05, 6.932032e-05, 7.320991e-05, 6.739325e-05,
                                                                                      3.838246e-05, 5.726806e-06, 7.094123e-07, 1.116956e-07],
                                                                                     [2.237953e-05, 5.685410e-05, 6.932032e-05, 7.320991e-05, 6.739325e-05,
                                                                                      3.838246e-05, 5.726806e-06, 7.094123e-07, 1.116956e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.409591e-05, 1.118391e-04, 1.361312e-04, 1.430246e-04, 1.306007e-04,
                                                                                      7.229734e-05, 1.013694e-05, 1.352289e-06, 2.105483e-07],
                                                                                     [4.409591e-05, 1.118391e-04, 1.361312e-04, 1.430246e-04, 1.306007e-04,
                                                                                      7.229734e-05, 1.013694e-05, 1.352289e-06, 2.105483e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.057208e-05, 1.030890e-04, 1.261760e-04, 1.337365e-04, 1.238674e-04,
                                                                                      7.227049e-05, 1.157070e-05, 1.321011e-06, 2.124075e-07],
                                                                                     [4.057208e-05, 1.030890e-04, 1.261760e-04, 1.337365e-04, 1.238674e-04,
                                                                                      7.227049e-05, 1.157070e-05, 1.321011e-06, 2.124075e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.017738e+00, 1.014514e+00, 1.010980e+00, 1.000546e+00, 9.844137e-01,
                                                                                   9.299548e-01, 8.212772e-01, 9.515456e-01, 9.310076e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.354480e-01, 9.357759e-01, 9.432074e-01, 9.499499e-01, 9.614999e-01,
                                                                                   1.009019e+00, 1.161774e+00, 9.883964e-01, 1.024863e+00],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.352452e-04, 5.006314e-04, 4.107292e-04, 3.891627e-04, 4.232981e-04,
                                                                                      4.963980e-04, 1.998284e-03, 5.385523e-03, 2.082378e-02],
                                                                                     [6.352452e-04, 5.006314e-04, 4.107292e-04, 3.891627e-04, 4.232981e-04,
                                                                                      4.963980e-04, 1.998284e-03, 5.385523e-03, 2.082378e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.251667e-03, 9.848043e-04, 8.065896e-04, 7.602774e-04, 8.203051e-04,
                                                                                      9.350172e-04, 3.537135e-03, 1.026594e-02, 3.925321e-02],
                                                                                     [1.251667e-03, 9.848043e-04, 8.065896e-04, 7.602774e-04, 8.203051e-04,
                                                                                      9.350172e-04, 3.537135e-03, 1.026594e-02, 3.925321e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.151643e-03, 9.077549e-04, 7.476042e-04, 7.109045e-04, 7.780132e-04,
                                                                                      9.346699e-04, 4.037424e-03, 1.002849e-02, 3.959983e-02],
                                                                                     [1.151643e-03, 9.077549e-04, 7.476042e-04, 7.109045e-04, 7.780132e-04,
                                                                                      9.346699e-04, 4.037424e-03, 1.002849e-02, 3.959983e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
