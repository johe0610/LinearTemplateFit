
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+02, 2.800000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02,
                                                                                   6.000000e+02, 6.800000e+02, 7.600000e+02, 8.400000e+02, 9.200000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+02, 2.800000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02,
                                                                                   6.000000e+02, 6.800000e+02, 7.600000e+02, 8.400000e+02, 9.200000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.000000e+02, 2.800000e+02, 3.600000e+02, 4.400000e+02, 5.200000e+02,
                                                                                   6.000000e+02, 6.800000e+02, 7.600000e+02, 8.400000e+02, 9.200000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.600000e+02, 2.400000e+02, 3.200000e+02, 4.000000e+02, 4.800000e+02,
                                                                                   5.600000e+02, 6.400000e+02, 7.200000e+02, 8.000000e+02, 8.800000e+02,
                                                                                   9.600000e+02],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.600000e+02, 2.400000e+02, 3.200000e+02, 4.000000e+02, 4.800000e+02,
                                                                                   5.600000e+02, 6.400000e+02, 7.200000e+02, 8.000000e+02, 8.800000e+02,
                                                                                   9.600000e+02],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.600000e+02, 2.400000e+02, 3.200000e+02, 4.000000e+02, 4.800000e+02,
                                                                                   5.600000e+02, 6.400000e+02, 7.200000e+02, 8.000000e+02, 8.800000e+02,
                                                                                   9.600000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.691946e-05, 2.894214e-04, 8.455826e-04, 9.552768e-04, 7.795680e-04,
                                                                                   5.721052e-04, 4.072150e-04, 2.882074e-04, 2.056282e-04, 1.492953e-04],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.756856e-05, 3.051098e-04, 8.728870e-04, 9.807763e-04, 7.966232e-04,
                                                                                   5.818144e-04, 4.124666e-04, 2.934442e-04, 2.087821e-04, 1.517406e-04],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.574513e-05, 2.767889e-04, 8.208010e-04, 9.397261e-04, 7.684676e-04,
                                                                                   5.638689e-04, 4.049532e-04, 2.878804e-04, 2.080948e-04, 1.467264e-04],
}

xerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                     [4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01,
                                                                                      4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01, 4.000000e+01],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.724665e-07, 7.131303e-07, 1.219036e-06, 1.296468e-06, 1.171915e-06,
                                                                                      1.004931e-06, 8.488382e-07, 7.149328e-07, 6.046843e-07, 5.154573e-07],
                                                                                     [1.724665e-07, 7.131303e-07, 1.219036e-06, 1.296468e-06, 1.171915e-06,
                                                                                      1.004931e-06, 8.488382e-07, 7.149328e-07, 6.046843e-07, 5.154573e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.776227e-07, 1.991069e-06, 3.369309e-06, 3.571248e-06, 3.222272e-06,
                                                                                      2.756598e-06, 2.321891e-06, 1.961594e-06, 1.656906e-06, 1.414857e-06],
                                                                                     [4.776227e-07, 1.991069e-06, 3.369309e-06, 3.571248e-06, 3.222272e-06,
                                                                                      2.756598e-06, 2.321891e-06, 1.961594e-06, 1.656906e-06, 1.414857e-06],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.340272e-07, 1.819985e-06, 3.135837e-06, 3.355788e-06, 3.036494e-06,
                                                                                      2.603214e-06, 2.207797e-06, 1.864004e-06, 1.586010e-06, 1.335139e-06],
                                                                                     [4.340272e-07, 1.819985e-06, 3.135837e-06, 3.355788e-06, 3.036494e-06,
                                                                                      2.603214e-06, 2.207797e-06, 1.864004e-06, 1.586010e-06, 1.335139e-06],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.038364e+00, 1.054206e+00, 1.032291e+00, 1.026693e+00, 1.021878e+00,
                                                                                   1.016971e+00, 1.012896e+00, 1.018170e+00, 1.015338e+00, 1.016379e+00],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.305929e-01, 9.563526e-01, 9.706929e-01, 9.837213e-01, 9.857608e-01,
                                                                                   9.856035e-01, 9.944457e-01, 9.988654e-01, 1.011995e+00, 9.827932e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410472.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.019338e-02, 2.463986e-03, 1.441652e-03, 1.357165e-03, 1.503288e-03,
                                                                                      1.756549e-03, 2.084496e-03, 2.480619e-03, 2.940668e-03, 3.452602e-03],
                                                                                     [1.019338e-02, 2.463986e-03, 1.441652e-03, 1.357165e-03, 1.503288e-03,
                                                                                      1.756549e-03, 2.084496e-03, 2.480619e-03, 2.940668e-03, 3.452602e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411053.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.822919e-02, 6.879481e-03, 3.984601e-03, 3.738443e-03, 4.133407e-03,
                                                                                      4.818341e-03, 5.701880e-03, 6.806189e-03, 8.057776e-03, 9.476902e-03],
                                                                                     [2.822919e-02, 6.879481e-03, 3.984601e-03, 3.738443e-03, 4.133407e-03,
                                                                                      4.818341e-03, 5.701880e-03, 6.806189e-03, 8.057776e-03, 9.476902e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411058.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.565254e-02, 6.288357e-03, 3.708493e-03, 3.512896e-03, 3.895098e-03,
                                                                                      4.550237e-03, 5.421699e-03, 6.467579e-03, 7.712999e-03, 8.942941e-03],
                                                                                     [2.565254e-02, 6.288357e-03, 3.708493e-03, 3.512896e-03, 3.895098e-03,
                                                                                      4.550237e-03, 5.421699e-03, 6.467579e-03, 7.712999e-03, 8.942941e-03],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
