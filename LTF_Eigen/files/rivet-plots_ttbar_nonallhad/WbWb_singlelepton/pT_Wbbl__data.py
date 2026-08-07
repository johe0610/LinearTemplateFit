
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 9.000000e+01, 1.500000e+02, 2.100000e+02, 2.700000e+02,
                                                                                   3.350000e+02, 4.100000e+02, 5.000000e+02, 6.250000e+02, 8.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 9.000000e+01, 1.500000e+02, 2.100000e+02, 2.700000e+02,
                                                                                   3.350000e+02, 4.100000e+02, 5.000000e+02, 6.250000e+02, 8.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [3.000000e+01, 9.000000e+01, 1.500000e+02, 2.100000e+02, 2.700000e+02,
                                                                                   3.350000e+02, 4.100000e+02, 5.000000e+02, 6.250000e+02, 8.000000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 6.000000e+01, 1.200000e+02, 1.800000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.700000e+02, 4.500000e+02, 5.500000e+02, 7.000000e+02,
                                                                                   9.000000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 6.000000e+01, 1.200000e+02, 1.800000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.700000e+02, 4.500000e+02, 5.500000e+02, 7.000000e+02,
                                                                                   9.000000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [0.000000e+00, 6.000000e+01, 1.200000e+02, 1.800000e+02, 2.400000e+02,
                                                                                   3.000000e+02, 3.700000e+02, 4.500000e+02, 5.500000e+02, 7.000000e+02,
                                                                                   9.000000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.171683e-01, 1.066176e-01, 3.417256e-02, 1.014558e-02, 3.460042e-03,
                                                                                   1.345664e-03, 5.729824e-04, 2.624680e-04, 1.133229e-04, 3.631189e-05],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.177390e-01, 1.053968e-01, 3.315825e-02, 9.756395e-03, 3.311277e-03,
                                                                                   1.281584e-03, 5.389750e-04, 2.517101e-04, 1.077643e-04, 3.501932e-05],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.128141e-01, 1.019928e-01, 3.248087e-02, 9.618343e-03, 3.271825e-03,
                                                                                   1.257979e-03, 5.281301e-04, 2.456691e-04, 1.094166e-04, 3.559618e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                     [3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01, 3.000000e+01,
                                                                                      3.500000e+01, 4.000000e+01, 5.000000e+01, 7.500000e+01, 1.000000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.331723e-05, 3.181052e-05, 1.805724e-05, 9.876236e-06, 5.789809e-06,
                                                                                      3.352498e-06, 2.047443e-06, 1.235607e-06, 6.594054e-07, 3.217927e-07],
                                                                                     [3.331723e-05, 3.181052e-05, 1.805724e-05, 9.876236e-06, 5.789809e-06,
                                                                                      3.352498e-06, 2.047443e-06, 1.235607e-06, 6.594054e-07, 3.217927e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.522848e-05, 6.177307e-05, 3.474370e-05, 1.891741e-05, 1.106431e-05,
                                                                                      6.399175e-06, 3.882640e-06, 2.361450e-06, 1.253752e-06, 6.171008e-07],
                                                                                     [6.522848e-05, 6.177307e-05, 3.474370e-05, 1.891741e-05, 1.106431e-05,
                                                                                      6.399175e-06, 3.882640e-06, 2.361450e-06, 1.253752e-06, 6.171008e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [6.127743e-05, 5.831845e-05, 3.299507e-05, 1.802055e-05, 1.056240e-05,
                                                                                      6.087131e-06, 3.689902e-06, 2.242946e-06, 1.214170e-06, 5.952535e-07],
                                                                                     [6.127743e-05, 5.831845e-05, 3.299507e-05, 1.802055e-05, 1.056240e-05,
                                                                                      6.087131e-06, 3.689902e-06, 2.242946e-06, 1.214170e-06, 5.952535e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.004871e+00, 9.885497e-01, 9.703180e-01, 9.616399e-01, 9.570049e-01,
                                                                                   9.523804e-01, 9.406484e-01, 9.590125e-01, 9.509490e-01, 9.644037e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.628381e-01, 9.566225e-01, 9.504957e-01, 9.480328e-01, 9.456027e-01,
                                                                                   9.348389e-01, 9.217213e-01, 9.359964e-01, 9.655295e-01, 9.802899e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [2.843536e-04, 2.983609e-04, 5.284134e-04, 9.734521e-04, 1.673335e-03,
                                                                                      2.491334e-03, 3.573309e-03, 4.707648e-03, 5.818819e-03, 8.861910e-03],
                                                                                     [2.843536e-04, 2.983609e-04, 5.284134e-04, 9.734521e-04, 1.673335e-03,
                                                                                      2.491334e-03, 3.573309e-03, 4.707648e-03, 5.818819e-03, 8.861910e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.567076e-04, 5.793891e-04, 1.016713e-03, 1.864596e-03, 3.197739e-03,
                                                                                      4.755403e-03, 6.776194e-03, 8.997097e-03, 1.106354e-02, 1.699446e-02],
                                                                                     [5.567076e-04, 5.793891e-04, 1.016713e-03, 1.864596e-03, 3.197739e-03,
                                                                                      4.755403e-03, 6.776194e-03, 8.997097e-03, 1.106354e-02, 1.699446e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [5.229864e-04, 5.469871e-04, 9.655428e-04, 1.776197e-03, 3.052680e-03,
                                                                                      4.523515e-03, 6.439817e-03, 8.545598e-03, 1.071425e-02, 1.639280e-02],
                                                                                     [5.229864e-04, 5.469871e-04, 9.655428e-04, 1.776197e-03, 3.052680e-03,
                                                                                      4.523515e-03, 6.439817e-03, 8.545598e-03, 1.071425e-02, 1.639280e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
