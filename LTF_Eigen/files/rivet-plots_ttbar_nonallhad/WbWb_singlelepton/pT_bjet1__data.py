
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda',
  'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda'
]

xpoints = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [4.250000e+01, 8.000000e+01, 1.200000e+02, 1.600000e+02, 2.050000e+02,
                                                                                   2.750000e+02, 3.800000e+02, 5.200000e+02, 7.900000e+02],
}

xedges = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [2.500000e+01, 6.000000e+01, 1.000000e+02, 1.400000e+02, 1.800000e+02,
                                                                                   2.300000e+02, 3.200000e+02, 4.400000e+02, 6.000000e+02, 9.800000e+02],
}

ref_xerrs = [
  [abs(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]   - xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))],
  [abs(xedges['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i+1] - xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda'][i]) for i in range(len(xpoints['user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda']))]
]

yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.297467e-02, 1.498704e-01, 9.406652e-02, 4.551008e-02, 1.923542e-02,
                                                                                   5.690591e-03, 1.127391e-03, 2.050680e-04, 2.986861e-05],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.523621e-02, 1.492602e-01, 9.236939e-02, 4.415280e-02, 1.842191e-02,
                                                                                   5.423772e-03, 1.066160e-03, 1.934652e-04, 2.895094e-05],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [8.766542e-02, 1.437192e-01, 9.053686e-02, 4.394357e-02, 1.856464e-02,
                                                                                   5.489490e-03, 1.083697e-03, 1.933307e-04, 2.739352e-05],
}

xerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                     [1.750000e+01, 2.000000e+01, 2.000000e+01, 2.000000e+01, 2.500000e+01,
                                                                                      4.500000e+01, 6.000000e+01, 8.000000e+01, 1.900000e+02],
                                                                                  ],
}

yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [3.881516e-05, 4.612731e-05, 3.660627e-05, 2.552617e-05, 1.489388e-05,
                                                                                      6.069371e-06, 2.360979e-06, 8.778840e-07, 2.165074e-07],
                                                                                     [3.881516e-05, 4.612731e-05, 3.660627e-05, 2.552617e-05, 1.489388e-05,
                                                                                      6.069371e-06, 2.360979e-06, 8.778840e-07, 2.165074e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.671526e-05, 8.991284e-05, 7.085098e-05, 4.911015e-05, 2.847961e-05,
                                                                                      1.157723e-05, 4.485445e-06, 1.665049e-06, 4.165947e-07],
                                                                                     [7.671526e-05, 8.991284e-05, 7.085098e-05, 4.911015e-05, 2.847961e-05,
                                                                                      1.157723e-05, 4.485445e-06, 1.665049e-06, 4.165947e-07],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.065375e-05, 8.466575e-05, 6.731049e-05, 4.700485e-05, 2.742269e-05,
                                                                                      1.117477e-05, 4.341442e-06, 1.600866e-06, 3.909912e-07],
                                                                                     [7.065375e-05, 8.466575e-05, 6.731049e-05, 4.700485e-05, 2.742269e-05,
                                                                                      1.117477e-05, 4.341442e-06, 1.600866e-06, 3.909912e-07],
                                                                                  ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                                                                   1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [1.024324e+00, 9.959285e-01, 9.819582e-01, 9.701763e-01, 9.577077e-01,
                                                                                   9.531123e-01, 9.456879e-01, 9.434197e-01, 9.692764e-01],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [9.428957e-01, 9.589565e-01, 9.624770e-01, 9.655788e-01, 9.651279e-01,
                                                                                   9.646608e-01, 9.612433e-01, 9.427639e-01, 9.171341e-01],
}

ratio0_yerrs = {
    'user.johessle.mc15_13TeV.410470.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [4.174810e-04, 3.077813e-04, 3.891530e-04, 5.608905e-04, 7.742945e-04,
                                                                                      1.066563e-03, 2.094197e-03, 4.280941e-03, 7.248660e-03],
                                                                                     [4.174810e-04, 3.077813e-04, 3.891530e-04, 5.608905e-04, 7.742945e-04,
                                                                                      1.066563e-03, 2.094197e-03, 4.280941e-03, 7.248660e-03],
                                                                                  ],
    'user.johessle.mc15_13TeV.411045.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [8.251200e-04, 5.999373e-04, 7.532008e-04, 1.079105e-03, 1.480582e-03,
                                                                                      2.034451e-03, 3.978606e-03, 8.119497e-03, 1.394758e-02],
                                                                                     [8.251200e-04, 5.999373e-04, 7.532008e-04, 1.079105e-03, 1.480582e-03,
                                                                                      2.034451e-03, 3.978606e-03, 8.119497e-03, 1.394758e-02],
                                                                                  ],
    'user.johessle.mc15_13TeV.411050.WbWb_singlelepton.test_topmass_mppui.yoda' : [
                                                                                     [7.599247e-04, 5.649264e-04, 7.155627e-04, 1.032845e-03, 1.425635e-03,
                                                                                      1.963727e-03, 3.850875e-03, 7.806513e-03, 1.309037e-02],
                                                                                     [7.599247e-04, 5.649264e-04, 7.155627e-04, 1.032845e-03, 1.425635e-03,
                                                                                      1.963727e-03, 3.850875e-03, 7.806513e-03, 1.309037e-02],
                                                                                  ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
