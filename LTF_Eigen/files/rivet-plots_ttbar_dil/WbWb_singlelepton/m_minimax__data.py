
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.700000e+02,
                                                                                   3.300000e+02, 4.000000e+02, 4.800000e+02, 5.700000e+02, 6.700000e+02,
                                                                                   8.100000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.700000e+02,
                                                                                   3.300000e+02, 4.000000e+02, 4.800000e+02, 5.700000e+02, 6.700000e+02,
                                                                                   8.100000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+02, 1.400000e+02, 1.800000e+02, 2.200000e+02, 2.700000e+02,
                                                                                   3.300000e+02, 4.000000e+02, 4.800000e+02, 5.700000e+02, 6.700000e+02,
                                                                                   8.100000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02, 6.200000e+02,
                                                                                   7.200000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02, 6.200000e+02,
                                                                                   7.200000e+02, 9.000000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.000000e+01, 1.200000e+02, 1.600000e+02, 2.000000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02, 6.200000e+02,
                                                                                   7.200000e+02, 9.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.837247e-04, 2.244575e-03, 2.311156e-03, 1.569037e-03, 9.017704e-04,
                                                                                   4.789438e-04, 2.552356e-04, 1.321645e-04, 6.944071e-05, 3.597604e-05,
                                                                                   1.592702e-05],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [6.095937e-04, 2.302194e-03, 2.376713e-03, 1.609142e-03, 9.154436e-04,
                                                                                   4.865487e-04, 2.590942e-04, 1.355695e-04, 7.199623e-05, 3.715326e-05,
                                                                                   1.554853e-05],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [5.658632e-04, 2.182575e-03, 2.274864e-03, 1.554064e-03, 8.924959e-04,
                                                                                   4.768313e-04, 2.548429e-04, 1.336533e-04, 6.895289e-05, 3.610834e-05,
                                                                                   1.571832e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                     [2.000000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 3.000000e+01,
                                                                                      3.000000e+01, 4.000000e+01, 4.000000e+01, 5.000000e+01, 5.000000e+01,
                                                                                      9.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.436083e-06, 2.814850e-06, 2.854813e-06, 2.352447e-06, 1.456665e-06,
                                                                                      1.062608e-06, 6.720130e-07, 4.838445e-07, 3.136198e-07, 2.259925e-07,
                                                                                      1.119251e-07],
                                                                                     [1.436083e-06, 2.814850e-06, 2.854813e-06, 2.352447e-06, 1.456665e-06,
                                                                                      1.062608e-06, 6.720130e-07, 4.838445e-07, 3.136198e-07, 2.259925e-07,
                                                                                      1.119251e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.990352e-06, 7.752777e-06, 7.873841e-06, 6.478258e-06, 3.990871e-06,
                                                                                      2.914228e-06, 1.842477e-06, 1.332672e-06, 8.694887e-07, 6.237917e-07,
                                                                                      3.008860e-07],
                                                                                     [3.990352e-06, 7.752777e-06, 7.873841e-06, 6.478258e-06, 3.990871e-06,
                                                                                      2.914228e-06, 1.842477e-06, 1.332672e-06, 8.694887e-07, 6.237917e-07,
                                                                                      3.008860e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.689276e-06, 7.244733e-06, 7.393336e-06, 6.107653e-06, 3.782392e-06,
                                                                                      2.765728e-06, 1.752449e-06, 1.270713e-06, 8.155211e-07, 5.908552e-07,
                                                                                      2.900245e-07],
                                                                                     [3.689276e-06, 7.244733e-06, 7.393336e-06, 6.107653e-06, 3.782392e-06,
                                                                                      2.765728e-06, 1.752449e-06, 1.270713e-06, 8.155211e-07, 5.908552e-07,
                                                                                      2.900245e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.044317e+00, 1.025670e+00, 1.028365e+00, 1.025560e+00, 1.015163e+00,
                                                                                   1.015878e+00, 1.015118e+00, 1.025763e+00, 1.036801e+00, 1.032722e+00,
                                                                                   9.762360e-01],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.694008e-01, 9.723778e-01, 9.842970e-01, 9.904572e-01, 9.897152e-01,
                                                                                   9.955893e-01, 9.984614e-01, 1.011265e+00, 9.929750e-01, 1.003677e+00,
                                                                                   9.868965e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.460206e-03, 1.254068e-03, 1.235232e-03, 1.499294e-03, 1.615339e-03,
                                                                                      2.218649e-03, 2.632912e-03, 3.660926e-03, 4.516368e-03, 6.281750e-03,
                                                                                      7.027372e-03],
                                                                                     [2.460206e-03, 1.254068e-03, 1.235232e-03, 1.499294e-03, 1.615339e-03,
                                                                                      2.218649e-03, 2.632912e-03, 3.660926e-03, 4.516368e-03, 6.281750e-03,
                                                                                      7.027372e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.836017e-03, 3.454007e-03, 3.406884e-03, 4.128811e-03, 4.425595e-03,
                                                                                      6.084697e-03, 7.218730e-03, 1.008343e-02, 1.252131e-02, 1.733909e-02,
                                                                                      1.889154e-02],
                                                                                     [6.836017e-03, 3.454007e-03, 3.406884e-03, 4.128811e-03, 4.425595e-03,
                                                                                      6.084697e-03, 7.218730e-03, 1.008343e-02, 1.252131e-02, 1.733909e-02,
                                                                                      1.889154e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.320233e-03, 3.227664e-03, 3.198977e-03, 3.892612e-03, 4.194407e-03,
                                                                                      5.774640e-03, 6.866005e-03, 9.614632e-03, 1.174414e-02, 1.642358e-02,
                                                                                      1.820959e-02],
                                                                                     [6.320233e-03, 3.227664e-03, 3.198977e-03, 3.892612e-03, 4.194407e-03,
                                                                                      5.774640e-03, 6.866005e-03, 9.614632e-03, 1.174414e-02, 1.642358e-02,
                                                                                      1.820959e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
