
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+02, 2.800000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02,
                                                                                   6.000000e+02, 6.800000e+02, 7.600000e+02, 8.400000e+02, 9.200000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+02, 2.800000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02,
                                                                                   6.000000e+02, 6.800000e+02, 7.600000e+02, 8.400000e+02, 9.200000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+02, 2.800000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02,
                                                                                   6.000000e+02, 6.800000e+02, 7.600000e+02, 8.400000e+02, 9.200000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.600000e+02, 2.400000e+02, 3.200000e+02, 4.000000e+02, 4.800000e+02,
                                                                                   5.600000e+02, 6.400000e+02, 7.200000e+02, 8.000000e+02, 8.800000e+02,
                                                                                   9.600000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.600000e+02, 2.400000e+02, 3.200000e+02, 4.000000e+02, 4.800000e+02,
                                                                                   5.600000e+02, 6.400000e+02, 7.200000e+02, 8.000000e+02, 8.800000e+02,
                                                                                   9.600000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.600000e+02, 2.400000e+02, 3.200000e+02, 4.000000e+02, 4.800000e+02,
                                                                                   5.600000e+02, 6.400000e+02, 7.200000e+02, 8.000000e+02, 8.800000e+02,
                                                                                   9.600000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.112485e-04, 1.978530e-02, 5.225540e-02, 4.625653e-02, 3.125676e-02,
                                                                                   1.979344e-02, 1.240354e-02, 7.853108e-03, 5.060566e-03, 3.326277e-03],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.364851e-04, 2.072701e-02, 5.257262e-02, 4.579136e-02, 3.064501e-02,
                                                                                   1.924953e-02, 1.204252e-02, 7.584056e-03, 4.862724e-03, 3.173759e-03],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.516819e-04, 1.817133e-02, 5.013072e-02, 4.477564e-02, 3.023845e-02,
                                                                                   1.905884e-02, 1.195191e-02, 7.518787e-03, 4.847808e-03, 3.159391e-03],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.483931e-06, 1.182960e-05, 1.923299e-05, 1.810941e-05, 1.490546e-05,
                                                                                      1.188352e-05, 9.433232e-06, 7.528061e-06, 6.063403e-06, 4.935463e-06],
                                                                                     [1.483931e-06, 1.182960e-05, 1.923299e-05, 1.810941e-05, 1.490546e-05,
                                                                                      1.188352e-05, 9.433232e-06, 7.528061e-06, 6.063403e-06, 4.935463e-06],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.011777e-06, 2.364651e-05, 3.767531e-05, 3.518919e-05, 2.882374e-05,
                                                                                      2.289252e-05, 1.815956e-05, 1.445315e-05, 1.162484e-05, 9.422923e-06],
                                                                                     [3.011777e-06, 2.364651e-05, 3.767531e-05, 3.518919e-05, 2.882374e-05,
                                                                                      2.289252e-05, 1.815956e-05, 1.445315e-05, 1.162484e-05, 9.422923e-06],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.501472e-06, 2.125062e-05, 3.530914e-05, 3.339441e-05, 2.747529e-05,
                                                                                      2.185649e-05, 1.735380e-05, 1.380559e-05, 1.112429e-05, 9.017397e-06],
                                                                                     [2.501472e-06, 2.125062e-05, 3.530914e-05, 3.339441e-05, 2.747529e-05,
                                                                                      2.185649e-05, 1.735380e-05, 1.380559e-05, 1.112429e-05, 9.017397e-06],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.081082e+00, 1.047596e+00, 1.006071e+00, 9.899437e-01, 9.804282e-01,
                                                                                   9.725207e-01, 9.708938e-01, 9.657394e-01, 9.609052e-01, 9.541475e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.086204e-01, 9.184258e-01, 9.593405e-01, 9.679853e-01, 9.674211e-01,
                                                                                   9.628867e-01, 9.635886e-01, 9.574282e-01, 9.579577e-01, 9.498280e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.767673e-03, 5.978984e-04, 3.680575e-04, 3.914995e-04, 4.768716e-04,
                                                                                      6.003767e-04, 7.605274e-04, 9.586091e-04, 1.198167e-03, 1.483780e-03],
                                                                                     [4.767673e-03, 5.978984e-04, 3.680575e-04, 3.914995e-04, 4.768716e-04,
                                                                                      6.003767e-04, 7.605274e-04, 9.586091e-04, 1.198167e-03, 1.483780e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [9.676439e-03, 1.195155e-03, 7.209841e-04, 7.607399e-04, 9.221602e-04,
                                                                                      1.156571e-03, 1.464063e-03, 1.840437e-03, 2.297142e-03, 2.832874e-03],
                                                                                     [9.676439e-03, 1.195155e-03, 7.209841e-04, 7.607399e-04, 9.221602e-04,
                                                                                      1.156571e-03, 1.464063e-03, 1.840437e-03, 2.297142e-03, 2.832874e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.036897e-03, 1.074061e-03, 6.757032e-04, 7.219394e-04, 8.790191e-04,
                                                                                      1.104229e-03, 1.399101e-03, 1.757978e-03, 2.198230e-03, 2.710958e-03],
                                                                                     [8.036897e-03, 1.074061e-03, 6.757032e-04, 7.219394e-04, 8.790191e-04,
                                                                                      1.104229e-03, 1.399101e-03, 1.757978e-03, 2.198230e-03, 2.710958e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
