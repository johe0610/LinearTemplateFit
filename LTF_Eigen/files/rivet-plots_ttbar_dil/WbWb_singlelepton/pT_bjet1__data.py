
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.441506e-03, 3.202151e-03, 2.513998e-03, 1.412674e-03, 6.692658e-04,
                                                                                   2.212590e-04, 5.122978e-05, 1.119286e-05, 1.458165e-06],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.500224e-03, 3.319727e-03, 2.571985e-03, 1.426845e-03, 6.704474e-04,
                                                                                   2.211722e-04, 5.143364e-05, 1.140624e-05, 1.485484e-06],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.377394e-03, 3.117294e-03, 2.489384e-03, 1.412577e-03, 6.738428e-04,
                                                                                   2.262499e-04, 5.152856e-05, 1.131967e-05, 1.453889e-06],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.407737e-06, 3.356904e-06, 2.977316e-06, 2.235628e-06, 1.378683e-06,
                                                                                      5.923997e-07, 2.477074e-07, 1.008984e-07, 2.368481e-08],
                                                                                     [2.407737e-06, 3.356904e-06, 2.977316e-06, 2.235628e-06, 1.378683e-06,
                                                                                      5.923997e-07, 2.477074e-07, 1.008984e-07, 2.368481e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.676714e-06, 9.296382e-06, 8.189627e-06, 6.108112e-06, 3.756424e-06,
                                                                                      1.612936e-06, 6.760545e-07, 2.766218e-07, 6.513069e-08],
                                                                                     [6.676714e-06, 9.296382e-06, 8.189627e-06, 6.108112e-06, 3.756424e-06,
                                                                                      1.612936e-06, 6.760545e-07, 2.766218e-07, 6.513069e-08],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.141375e-06, 8.645511e-06, 7.731258e-06, 5.832402e-06, 3.608477e-06,
                                                                                      1.562551e-06, 6.502136e-07, 2.653813e-07, 6.186621e-08],
                                                                                     [6.141375e-06, 8.645511e-06, 7.731258e-06, 5.832402e-06, 3.608477e-06,
                                                                                      1.562551e-06, 6.502136e-07, 2.653813e-07, 6.186621e-08],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.040734e+00, 1.036718e+00, 1.023066e+00, 1.010031e+00, 1.001766e+00,
                                                                                   9.996077e-01, 1.003979e+00, 1.019064e+00, 1.018735e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.555243e-01, 9.735000e-01, 9.902092e-01, 9.999313e-01, 1.006839e+00,
                                                                                   1.022557e+00, 1.005832e+00, 1.011330e+00, 9.970675e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.670293e-03, 1.048328e-03, 1.184295e-03, 1.582551e-03, 2.059993e-03,
                                                                                      2.677404e-03, 4.835223e-03, 9.014532e-03, 1.624289e-02],
                                                                                     [1.670293e-03, 1.048328e-03, 1.184295e-03, 1.582551e-03, 2.059993e-03,
                                                                                      2.677404e-03, 4.835223e-03, 9.014532e-03, 1.624289e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.631763e-03, 2.903168e-03, 3.257611e-03, 4.323794e-03, 5.612754e-03,
                                                                                      7.289810e-03, 1.319651e-02, 2.471413e-02, 4.466620e-02],
                                                                                     [4.631763e-03, 2.903168e-03, 3.257611e-03, 4.323794e-03, 5.612754e-03,
                                                                                      7.289810e-03, 1.319651e-02, 2.471413e-02, 4.466620e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.260388e-03, 2.699907e-03, 3.075284e-03, 4.128626e-03, 5.391695e-03,
                                                                                      7.062090e-03, 1.269210e-02, 2.370987e-02, 4.242744e-02],
                                                                                     [4.260388e-03, 2.699907e-03, 3.075284e-03, 4.128626e-03, 5.391695e-03,
                                                                                      7.062090e-03, 1.269210e-02, 2.370987e-02, 4.242744e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
