
import numpy as np
from numpy import nan

add_legend_handle = [
    'WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz',
    'WbWb_LO_PS_Had_162/Analysis/WbWb_Slurm_Template_162_5.yoda.gz',
    'WbWb_LO_PS_Had_165/Analysis/WbWb_Slurm_Template_165.yoda.gz',
    'WbWb_LO_PS_Had_167/Analysis/WbWb_Slurm_Template_167_5.yoda.gz',
    'WbWb_LO_PS_Had_170/Analysis/WbWb_Slurm_Template_170.yoda.gz',
    'WbWb_LO_PS_Had_172/Analysis/WbWb_Slurm_Template_172_5.yoda.gz',
    'WbWb_LO_PS_Had_175/Analysis/WbWb_Slurm_Template_175.yoda.gz',
    'WbWb_LO_PS_Had_177/Analysis/WbWb_Slurm_Template_177_5.yoda.gz',
    'WbWb_LO_PS_Had_180/Analysis/WbWb_Slurm_Template_180.yoda.gz',
    'WbWb_LO_PS_Had_182/Analysis/WbWb_Slurm_Template_182_5.yoda.gz'
]

xpoints = {
    'WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_162/Analysis/WbWb_Slurm_Template_162_5.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_165/Analysis/WbWb_Slurm_Template_165.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_167/Analysis/WbWb_Slurm_Template_167_5.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_170/Analysis/WbWb_Slurm_Template_170.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_172/Analysis/WbWb_Slurm_Template_172_5.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_175/Analysis/WbWb_Slurm_Template_175.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_177/Analysis/WbWb_Slurm_Template_177_5.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_180/Analysis/WbWb_Slurm_Template_180.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
    'WbWb_LO_PS_Had_182/Analysis/WbWb_Slurm_Template_182_5.yoda.gz' : [40.0, 70.0, 102.5, 140.0, 380.0],
}
xedges = {
    'WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_162/Analysis/WbWb_Slurm_Template_162_5.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_165/Analysis/WbWb_Slurm_Template_165.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_167/Analysis/WbWb_Slurm_Template_167_5.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_170/Analysis/WbWb_Slurm_Template_170.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_172/Analysis/WbWb_Slurm_Template_172_5.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_175/Analysis/WbWb_Slurm_Template_175.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_177/Analysis/WbWb_Slurm_Template_177_5.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_180/Analysis/WbWb_Slurm_Template_180.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
    'WbWb_LO_PS_Had_182/Analysis/WbWb_Slurm_Template_182_5.yoda.gz' : [25.0, 55.0, 85.0, 120.0, 160.0, 600.0],
}
ref_xerrs = [
  [abs(xpoints['WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz'][i]   - xedges['WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz'][i]) for i in range(len(xpoints['WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz']))],
  [abs(xedges['WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz'][i+1] - xpoints['WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz'][i]) for i in range(len(xpoints['WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz']))]
]

yvals = {
    'WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz' : [0.04427853, 0.003626336, 0.0007609218, 0.0002821365, 2.275569e-05],
    'WbWb_LO_PS_Had_162/Analysis/WbWb_Slurm_Template_162_5.yoda.gz' : [0.04874602, 0.003767786, 0.0008003, 0.0002785069, 2.3427e-05],
    'WbWb_LO_PS_Had_165/Analysis/WbWb_Slurm_Template_165.yoda.gz' : [0.04760965, 0.003762608, 0.0007715767, 0.0002829549, 2.091129e-05],
    'WbWb_LO_PS_Had_167/Analysis/WbWb_Slurm_Template_167_5.yoda.gz' : [0.04635853, 0.00362606, 0.0008222653, 0.0002422707, 2.310784e-05],
    'WbWb_LO_PS_Had_170/Analysis/WbWb_Slurm_Template_170.yoda.gz' : [0.04511767, 0.003678725, 0.0007567258, 0.0002909821, 2.260942e-05],
    'WbWb_LO_PS_Had_172/Analysis/WbWb_Slurm_Template_172_5.yoda.gz' : [0.04440217, 0.003642909, 0.0007194706, 0.0002844227, 2.204094e-05],
    'WbWb_LO_PS_Had_175/Analysis/WbWb_Slurm_Template_175.yoda.gz' : [0.04281659, 0.003554709, 0.0007162226, 0.0002700909, 2.243363e-05],
    'WbWb_LO_PS_Had_177/Analysis/WbWb_Slurm_Template_177_5.yoda.gz' : [0.04126507, 0.003528139, 0.0007675302, 0.0002533275, 2.456048e-05],
    'WbWb_LO_PS_Had_180/Analysis/WbWb_Slurm_Template_180.yoda.gz' : [0.04038591, 0.003432569, 0.000747575, 0.0003009276, 2.315671e-05],
    'WbWb_LO_PS_Had_182/Analysis/WbWb_Slurm_Template_182_5.yoda.gz' : [0.03910919, 0.003408806, 0.0006977466, 0.0002629419, 2.213683e-05],
}
xerrs = {
    'WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_162/Analysis/WbWb_Slurm_Template_162_5.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_165/Analysis/WbWb_Slurm_Template_165.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_167/Analysis/WbWb_Slurm_Template_167_5.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_170/Analysis/WbWb_Slurm_Template_170.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_172/Analysis/WbWb_Slurm_Template_172_5.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_175/Analysis/WbWb_Slurm_Template_175.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_177/Analysis/WbWb_Slurm_Template_177_5.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_180/Analysis/WbWb_Slurm_Template_180.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
    'WbWb_LO_PS_Had_182/Analysis/WbWb_Slurm_Template_182_5.yoda.gz' : [
        [15.0, 15.0, 17.5, 20.0, 220.0],
        [15.0, 15.0, 17.5, 20.0, 220.0],
    ],
}
yerrs = {
    'WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz' : [
        [0.0001634456, 4.585937e-05, 1.927821e-05, 1.079205e-05, 9.134667e-07],
        [0.0001634456, 4.585937e-05, 1.927821e-05, 1.079205e-05, 9.134667e-07],
    ],
    'WbWb_LO_PS_Had_162/Analysis/WbWb_Slurm_Template_162_5.yoda.gz' : [
        [0.0001941324, 5.396314e-05, 2.937754e-05, 1.428858e-05, 1.178672e-06],
        [0.0001941324, 5.396314e-05, 2.937754e-05, 1.428858e-05, 1.178672e-06],
    ],
    'WbWb_LO_PS_Had_165/Analysis/WbWb_Slurm_Template_165.yoda.gz' : [
        [0.0001830256, 5.080283e-05, 2.060149e-05, 1.140977e-05, 9.109992e-07],
        [0.0001830256, 5.080283e-05, 2.060149e-05, 1.140977e-05, 9.109992e-07],
    ],
    'WbWb_LO_PS_Had_167/Analysis/WbWb_Slurm_Template_167_5.yoda.gz' : [
        [0.000372269, 0.0001014393, 4.529416e-05, 1.981211e-05, 2.055088e-06],
        [0.000372269, 0.0001014393, 4.529416e-05, 1.981211e-05, 2.055088e-06],
    ],
    'WbWb_LO_PS_Had_170/Analysis/WbWb_Slurm_Template_170.yoda.gz' : [
        [0.0001656729, 4.697977e-05, 1.890214e-05, 1.139776e-05, 9.0715e-07],
        [0.0001656729, 4.697977e-05, 1.890214e-05, 1.139776e-05, 9.0715e-07],
    ],
    'WbWb_LO_PS_Had_172/Analysis/WbWb_Slurm_Template_172_5.yoda.gz' : [
        [0.000189614, 5.105353e-05, 1.952383e-05, 1.095933e-05, 8.612971e-07],
        [0.000189614, 5.105353e-05, 1.952383e-05, 1.095933e-05, 8.612971e-07],
    ],
    'WbWb_LO_PS_Had_175/Analysis/WbWb_Slurm_Template_175.yoda.gz' : [
        [0.0001490132, 4.277409e-05, 2.042809e-05, 9.959789e-06, 9.696813e-07],
        [0.0001490132, 4.277409e-05, 2.042809e-05, 9.959789e-06, 9.696813e-07],
    ],
    'WbWb_LO_PS_Had_177/Analysis/WbWb_Slurm_Template_177_5.yoda.gz' : [
        [0.0002242446, 6.445984e-05, 2.856611e-05, 1.411852e-05, 1.445473e-06],
        [0.0002242446, 6.445984e-05, 2.856611e-05, 1.411852e-05, 1.445473e-06],
    ],
    'WbWb_LO_PS_Had_180/Analysis/WbWb_Slurm_Template_180.yoda.gz' : [
        [0.0001553709, 4.583196e-05, 1.864972e-05, 1.166222e-05, 9.180234e-07],
        [0.0001553709, 4.583196e-05, 1.864972e-05, 1.166222e-05, 9.180234e-07],
    ],
    'WbWb_LO_PS_Had_182/Analysis/WbWb_Slurm_Template_182_5.yoda.gz' : [
        [0.0001274432, 3.735764e-05, 1.55876e-05, 8.591634e-06, 7.71126e-07],
        [0.0001274432, 3.735764e-05, 1.55876e-05, 8.591634e-06, 7.71126e-07],
    ],
}
variation_yvals = {
}


# lists for ratio plot
ratio0_yvals = {
    'WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz' : [1.0, 1.0, 1.0, 1.0, 1.0],
    'WbWb_LO_PS_Had_162/Analysis/WbWb_Slurm_Template_162_5.yoda.gz' : [1.100895174252623, 1.039006313810965, 1.0517506529580307, 0.9871353050739624, 1.0295007534379315],
    'WbWb_LO_PS_Had_165/Analysis/WbWb_Slurm_Template_165.yoda.gz' : [1.0752310431263188, 1.0375784262682775, 1.0140026215571691, 1.0029007235859237, 0.9189477444981893],
    'WbWb_LO_PS_Had_167/Analysis/WbWb_Slurm_Template_167_5.yoda.gz' : [1.0469753625515572, 0.9999238901193932, 1.0806173512179569, 0.8587003099563508, 1.0154752503659523],
    'WbWb_LO_PS_Had_170/Analysis/WbWb_Slurm_Template_170.yoda.gz' : [1.0189513969863044, 1.0144468135329985, 0.9944856357118432, 1.0313522000875461, 0.9935721571176264],
    'WbWb_LO_PS_Had_172/Analysis/WbWb_Slurm_Template_172_5.yoda.gz' : [1.0027923239547472, 1.0045701777220866, 0.9455250197852132, 1.0081031699195249, 0.9685902734656695],
    'WbWb_LO_PS_Had_175/Analysis/WbWb_Slurm_Template_175.yoda.gz' : [0.966983095418931, 0.9802481071803606, 0.9412565128243139, 0.9573057722060067, 0.9858470562747164],
    'WbWb_LO_PS_Had_177/Analysis/WbWb_Slurm_Template_177_5.yoda.gz' : [0.9319430884448964, 0.9729211523697749, 1.008684729495199, 0.897889851189052, 1.0793115919578795],
    'WbWb_LO_PS_Had_180/Analysis/WbWb_Slurm_Template_180.yoda.gz' : [0.9120878674156525, 0.9465667274074989, 0.9824596955955264, 1.0666028677608177, 1.0176228450994016],
    'WbWb_LO_PS_Had_182/Analysis/WbWb_Slurm_Template_182_5.yoda.gz' : [0.8832540285325642, 0.940013832143519, 0.9169754368977208, 0.9319669734330723, 0.9728041645847698],
}
ratio0_yerrs = {
    'WbWb_LO_PS_Had_171/Analysis/WbWb_Slurm_Template_171.yoda.gz' : [
        [0.003691305921854226, 0.012646199910874227, 0.02533533669294269, 0.03825116565917561, 0.040142342420730816],
        [0.003691305921854226, 0.012646199910874227, 0.02533533669294269, 0.03825116565917561, 0.040142342420730816],
    ],
    'WbWb_LO_PS_Had_162/Analysis/WbWb_Slurm_Template_162_5.yoda.gz' : [
        [0.004384346092790343, 0.014880899067267899, 0.03860783066012828, 0.0506442094518079, 0.05179680334896459],
        [0.004384346092790343, 0.014880899067267899, 0.03860783066012828, 0.0506442094518079, 0.05179680334896459],
    ],
    'WbWb_LO_PS_Had_165/Analysis/WbWb_Slurm_Template_165.yoda.gz' : [
        [0.0041335066904885955, 0.014009410600672415, 0.027074385304771135, 0.04044060233255888, 0.04003390800278964],
        [0.0041335066904885955, 0.014009410600672415, 0.027074385304771135, 0.04044060233255888, 0.04003390800278964],
    ],
    'WbWb_LO_PS_Had_167/Analysis/WbWb_Slurm_Template_167_5.yoda.gz' : [
        [0.00840743809697386, 0.027972945695048667, 0.059525380926134595, 0.07022171891974274, 0.09031095079955827],
        [0.00840743809697386, 0.027972945695048667, 0.059525380926134595, 0.07022171891974274, 0.09031095079955827],
    ],
    'WbWb_LO_PS_Had_170/Analysis/WbWb_Slurm_Template_170.yoda.gz' : [
        [0.003741607953109554, 0.012955161904467759, 0.024841107193932413, 0.04039803428482312, 0.039864754705306676],
        [0.003741607953109554, 0.012955161904467759, 0.024841107193932413, 0.04039803428482312, 0.039864754705306676],
    ],
    'WbWb_LO_PS_Had_172/Analysis/WbWb_Slurm_Template_172_5.yoda.gz' : [
        [0.004282301151370653, 0.014078543742223556, 0.02565812938990577, 0.038844070157530136, 0.03784974659085266],
        [0.004282301151370653, 0.014078543742223556, 0.02565812938990577, 0.038844070157530136, 0.03784974659085266],
    ],
    'WbWb_LO_PS_Had_175/Analysis/WbWb_Slurm_Template_175.yoda.gz' : [
        [0.003365360141811392, 0.011795401749865429, 0.02684650380630441, 0.03530131337136457, 0.04261269598944264],
        [0.003365360141811392, 0.011795401749865429, 0.02684650380630441, 0.03530131337136457, 0.04261269598944264],
    ],
    'WbWb_LO_PS_Had_177/Analysis/WbWb_Slurm_Template_177_5.yoda.gz' : [
        [0.005064409319821593, 0.01777547364612656, 0.037541453011334415, 0.0500414515668834, 0.06352138739805298],
        [0.005064409319821593, 0.01777547364612656, 0.037541453011334415, 0.0500414515668834, 0.06352138739805298],
    ],
    'WbWb_LO_PS_Had_180/Analysis/WbWb_Slurm_Template_180.yoda.gz' : [
        [0.0035089444026258324, 0.012638641317296577, 0.0245093779676177, 0.04133538198708781, 0.04034258684311485],
        [0.0035089444026258324, 0.012638641317296577, 0.0245093779676177, 0.04133538198708781, 0.04034258684311485],
    ],
    'WbWb_LO_PS_Had_182/Analysis/WbWb_Slurm_Template_182_5.yoda.gz' : [
        [0.002878216598428177, 0.010301759130979589, 0.020485153664936397, 0.030452047147391425, 0.03388717283457456],
        [0.002878216598428177, 0.010301759130979589, 0.020485153664936397, 0.030452047147391425, 0.03388717283457456],
    ],
}
ratio0_variation_vals = {
}
ratio_band_edges = {
}