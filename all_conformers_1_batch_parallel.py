import rdkit 
from rdkit import Chem
from rdkit.Chem import AllChem
# from rdkit.Chem import SaltRemover
import glob
import cPickle
import numpy as np
from timeout import timeout
import time
import cProfile
from SmallestEnclosingCircle_CASTING import returnCircleAsTuple
import scipy
from sklearn.decomposition import PCA
import math, random, cmath
import sys
import multiprocessing
from subprocess import call
from miniball_example_containers import doit
import logging

""" Before running
- Build the returnCircleAsTuple library and module using pybind11
- tar the result
- Change origin of the original SDFile
- Launch the job with multiple cores, maybe two?
- split adequately the workload
"""

"""S=NC(=O)C(=O)N=S forms an exclusion complex at 7.9 A of miniball radius"""

n_steps = 100000
tol = 1E-9
n_threads = 1

# fname = '/home/macenrola/Thesis/exclusion_complexes/weakly_binding_pubchem'
# fname = '/home/macenrola/Thesis/test_rolling_on_server/Compound_000000001_000025000.sdf.gz.smi.can'
flist  = sys.argv[1]
flist = flist.split(':')
print flist
# print '{} is the proc number'.format(multiprocessing.cpu_count())

# def strip_salts(smi):
#     """Removes the salts from a molecule by removing the fragments preceded by a '.'"""
#     smi = smi.replace('\\', '\\\\')
#     salts = smi.split('.')
#     smi_no_salt = ''
#     for salt in salts:
#         if len(salt)>len(smi_no_salt):
#             smi_no_salt = salt
#     return smi_no_salt

def get_CB_BLOCK():
    """
    PRE: -
    POST: Returns the RDKIT_BLOCK of CB7 as a string
    """
    return 'CB7\n     RDKit          3D\n\n126147  0  0  0  0  0  0  0  0999 V2000\n    1.1730   -4.4176    2.9511 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0750   -4.6349    3.6688 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0765   -3.4911    4.7539 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1723   -2.7865    4.5031 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9191   -3.3583    3.4714 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3124   -4.4064    2.9369 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0522   -3.3353    3.4421 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3130   -2.7704    4.4841 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0691   -3.0468    3.1437 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1939   -3.0104    3.0984 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1854   -0.4434    5.2692 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0591   -0.0113    5.8877 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0441    1.5519    5.6863 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2039    1.7887    4.9742 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9401    0.6201    4.7698 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2994   -0.4233    5.2450 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0305    0.6539    4.7401 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2819    1.8110    4.9656 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0894    0.5504    4.3212 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1741    0.6051    4.2737 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2527   -3.8136   -3.7077 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0124   -3.8818   -4.4671 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0312   -2.5540   -5.3160 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2779   -1.9153   -4.9189 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0091   -2.6783   -4.0054 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2331   -3.7867   -3.7189 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9639   -2.6375   -4.0281 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2071   -1.8908   -4.9333 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1545   -2.4399   -3.6082 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1102   -2.3785   -3.6467 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3096    4.4630   -2.8373 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0731    5.2303   -2.7776 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0598    5.7883   -1.3030 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2961    5.2698   -0.7356 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0425    4.5245   -1.6504 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1763    4.4888   -2.8686 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9303    4.5503   -1.6948 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1898    5.2750   -0.7590 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1845    4.0859   -1.4747 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0830    4.1303   -1.5454 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2378    3.8620    3.6411 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0016    4.6239    3.6943 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0176    5.4457    2.3491 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2631    5.0331    1.7185 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9935    4.1349    2.4992 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2476    3.8752    3.6162 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9787    4.1602    2.4607 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2224    5.0509    1.6961 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1395    3.7339    2.2694 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1253    3.7745    2.2088 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8611   -0.7617   -5.5808 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8794   -1.7513    5.3501 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7603   -0.7297   -5.6064 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7412   -1.7795    5.3807 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1982   -5.2981    0.6479 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0498   -5.9266    0.2399 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0334   -5.7852   -1.3302 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2138   -5.0813   -1.5932 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9522   -4.8379   -0.4333 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2882   -5.2657    0.6262 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0214   -4.8068   -0.4705 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2708   -5.0748   -1.6171 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1034   -4.3917   -0.3828 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.1644   -4.3384   -0.4391 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7455   -5.3457    1.9916 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8706   -5.3200    1.9558 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8157   -4.8707   -2.9471 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7943   -4.9008   -2.9121 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7863    3.0997    4.7479 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8328    3.1315    4.7180 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8453    5.6668    0.5485 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7668    5.6860    0.5091 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1731    2.6985   -4.5639 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0766    2.5733   -5.2995 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0673    1.0736   -5.7843 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1835    0.5544   -5.2506 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.9264    1.5226   -4.5717 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3114    2.6728   -4.5336 N   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0441    1.4838   -4.5307 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3002    0.5315   -5.2313 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0788    1.3907   -4.1449 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1843    1.3310   -4.0791 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7319    3.9561   -4.0995 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8841    3.9189   -4.0557 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0817   -5.6424    4.1014 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0872   -3.8701    5.7828 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0702   -0.3130    6.9419 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0393    2.1106    6.6299 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0054   -4.7893   -5.0823 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0345   -2.7314   -6.3981 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0901    6.0195   -3.5388 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0617    6.8837   -1.2538 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0074    5.2580    4.5891 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0242    6.5312    2.5046 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7799   -0.8998   -6.6666 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9121   -0.7356   -5.2918 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9301   -1.6557    5.0751 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7985   -2.0876    6.3912 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8160   -0.6826   -5.3367 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6624   -0.8715   -6.6906 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.7974   -1.7001    5.1217 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6395   -2.1199    6.4192 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0657   -6.9691    0.5804 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0296   -6.7488   -1.8530 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8037   -5.0954    1.9086 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6356   -6.3657    2.3808 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9239   -5.0599    1.8472 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7805   -6.3449    2.3374 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8721   -4.6322   -2.8198 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7133   -5.8028   -3.5171 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8522   -4.6814   -2.7612 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6872   -5.8367   -3.4748 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8415    2.9390    4.5246 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6891    3.6916    5.6670 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8862    2.9907    4.4754 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7409    3.7289    5.6339 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9003    5.3920    0.5370 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7470    6.7545    0.6523 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6574    6.7728    0.6099 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8247    5.4251    0.4716 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0949    3.2963   -6.1235 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0742    0.9704   -6.8759 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7919    3.7784   -3.9142 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6170    4.7060   -4.8933 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8027    4.6722   -4.8487 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9349    3.7196   -3.8436 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1 65  1  0\n  2  3  1  0\n  2  6  1  0\n  2 85  1  0\n  3  4  1  0\n  3 86  1  0\n  4  5  1  0\n  4 54  1  0\n  5  1  1  0\n  6  7  1  0\n  6 66  1  0\n  7  8  1  0\n  8  3  1  0\n  9  5  2  0\n 10  7  2  0\n 11 12  1  0\n 12 13  1  0\n 12 16  1  0\n 12 87  1  0\n 13 14  1  0\n 13 88  1  0\n 14 15  1  0\n 14 69  1  0\n 15 11  1  0\n 16 17  1  0\n 17 18  1  0\n 18 13  1  0\n 18 70  1  0\n 19 15  2  0\n 20 17  2  0\n 21 22  1  0\n 21 68  1  0\n 22 23  1  0\n 22 26  1  0\n 22 89  1  0\n 23 24  1  0\n 23 90  1  0\n 24 25  1  0\n 25 21  1  0\n 26 27  1  0\n 26 67  1  0\n 27 28  1  0\n 28 23  1  0\n 28 53  1  0\n 29 25  2  0\n 30 27  2  0\n 31 32  1  0\n 31 84  1  0\n 32 33  1  0\n 32 36  1  0\n 32 91  1  0\n 33 34  1  0\n 33 92  1  0\n 34 35  1  0\n 34 71  1  0\n 35 31  1  0\n 36 37  1  0\n 36 83  1  0\n 37 38  1  0\n 38 33  1  0\n 38 72  1  0\n 39 35  2  0\n 40 37  2  0\n 41 42  1  0\n 41 69  1  0\n 42 43  1  0\n 42 46  1  0\n 42 93  1  0\n 43 44  1  0\n 43 94  1  0\n 44 45  1  0\n 45 41  1  0\n 46 47  1  0\n 46 70  1  0\n 47 48  1  0\n 48 43  1  0\n 48 72  1  0\n 49 45  2  0\n 50 47  2  0\n 51 24  1  0\n 51 80  1  0\n 51 95  1  0\n 51 96  1  0\n 52  8  1  0\n 52 16  1  0\n 52 97  1  0\n 52 98  1  0\n 53 99  1  0\n 53100  1  0\n 54 11  1  0\n 54101  1  0\n 54102  1  0\n 55 56  1  0\n 55 65  1  0\n 56 57  1  0\n 56 60  1  0\n 56103  1  0\n 57 58  1  0\n 57104  1  0\n 58 59  1  0\n 58 68  1  0\n 59 55  1  0\n 60 61  1  0\n 60 66  1  0\n 61 62  1  0\n 62 57  1  0\n 62 67  1  0\n 63 59  2  0\n 64 61  2  0\n 65105  1  0\n 65106  1  0\n 66107  1  0\n 66108  1  0\n 67109  1  0\n 67110  1  0\n 68111  1  0\n 68112  1  0\n 69113  1  0\n 69114  1  0\n 70115  1  0\n 70116  1  0\n 71 44  1  0\n 71117  1  0\n 71118  1  0\n 72119  1  0\n 72120  1  0\n 73 74  1  0\n 74 75  1  0\n 74 78  1  0\n 74121  1  0\n 75 76  1  0\n 75122  1  0\n 76 53  1  0\n 76 77  1  0\n 77 73  1  0\n 77 81  2  0\n 78 79  1  0\n 79 80  1  0\n 80 75  1  0\n 82 79  2  0\n 83 73  1  0\n 83123  1  0\n 83124  1  0\n 84 78  1  0\n 84125  1  0\n 84126  1  0\nM  END\n'


def get_atoms_coords(RDKIT_BLOCK):
    """Takes as input an RDKIT BLOCK and returns a list of atoms with a numpy array containing the coordinates"""
    RDKIT_BLOCK = RDKIT_BLOCK.split('\n')
    atm_number = int(RDKIT_BLOCK[3][:3])
    RDKIT_BLOCK = [x.split() for x in RDKIT_BLOCK]
    atm_list = []
    coords_array = np.zeros([atm_number, 3], dtype=float)
    for i, line in enumerate(RDKIT_BLOCK[4:4+atm_number]):
        coords_atm = line
        atm_list.append(coords_atm[3])
        coords_array[i, :] = coords_atm[:3]
    return atm_list, coords_array

# def generate_mol_from_MDL(mol, coord_matrix, infile):
#     """Will write the MDL of the mol file then replace the xyz coordinates from the coord_matrix"""
#     w = Chem.SDWriter(infile+'_tmp_mol.sdf')
#     w.write(mol)
#     w.close()
#     with open(infile+'_tmp_mol.sdf', 'r') as r:
#         RDKIT_BLOCK = r.readlines()
#         atm_number = int(RDKIT_BLOCK[3][:3])
#         for i in range(0,atm_number):
#             j = i+4
#             RDKIT_BLOCK[j] = RDKIT_BLOCK[j].split()
#             RDKIT_BLOCK[j][:3] = coord_matrix[i, :]
#             RDKIT_BLOCK[j] = (' '*(3+int(np.sign(RDKIT_BLOCK[j][0])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][0])+
#                             ' '*(3+int(np.sign(RDKIT_BLOCK[j][1])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][1])+
#                             ' '*(3+int(np.sign(RDKIT_BLOCK[j][2])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][2])+
#                             ' {}   '.format(RDKIT_BLOCK[j][3]) + '  '.join(RDKIT_BLOCK[j][4:]) + '\n'
#                             )

#         with open(infile+'_tmp_mol_align.sdf', 'wb') as w: 
#             w.write(''.join(RDKIT_BLOCK))
#     supp = Chem.SDMolSupplier(infile+'_tmp_mol_align.sdf', removeHs = False)
#     mol_align = supp[0]
#     return mol_align  

def generate_mol_from_MDL(MOL_IN, coord_matrix):
    """Will write the MDL of the mol file then replace the xyz coordinates from the coord_matrix"""
    RDKIT_BLOCK_IN = Chem.MolToMolBlock(MOL_IN)
    RDKIT_BLOCK = [x+'\n' for x in RDKIT_BLOCK_IN.split('\n')]
    atm_number = int(RDKIT_BLOCK[3][:3])
    for i in range(0,atm_number):
        j = i+4
        RDKIT_BLOCK[j] = RDKIT_BLOCK[j].split()
        RDKIT_BLOCK[j][:3] = coord_matrix[i, :]
        RDKIT_BLOCK[j] = (' '*(3+int(np.sign(RDKIT_BLOCK[j][0])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][0])+
                ' '*(3+int(np.sign(RDKIT_BLOCK[j][1])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][1])+
                ' '*(3+int(np.sign(RDKIT_BLOCK[j][2])==1)) + '{0:.4f}'.format(RDKIT_BLOCK[j][2])+
                ' {}   '.format(RDKIT_BLOCK[j][3]) + '  '.join(RDKIT_BLOCK[j][4:]) + '\n'
                )

    RDKIT_BLOCK_OUT = ''.join(RDKIT_BLOCK)

    return Chem.MolFromMolBlock(RDKIT_BLOCK_OUT, removeHs=False)  

# def get_volume(mol, target_mol):
#     """
#     Returns the volume of a molecule described by a SMILES as defined in the article
#     Zhao, Yuan H., Michael H. Abraham, and Andreas M. Zissimos. "Fast calculation of van der Waals volume as a sum of atomic and bond contributions and its application to drug compounds." The Journal of organic chemistry 68.19 (2003): 7368-7373.
#     REQUIRES: 
#     - OPENBABEL LIBRARIES
#     - Python Bindings
#     - string module
#     """
#     volume = 0
#     index_of_vols = {1:7.24, 6:20.58, 7:15.60, 8:14.71, 9:13.31, 17:22.45, 35:26.52, 53:32.52, 15:24.43, 16:24.43, 33:26.52, 5:40.48, 14:38.79, 34:28.73, 52:36.62} 
#     try:
#         for atom in mol.GetAtoms():
#             volume+=index_of_vols[atom.GetAtomicNum()]
#     except KeyError:
#         return 0
#     r_a  = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
#     r_na = Chem.rdMolDescriptors.CalcNumAliphaticRings(mol)
#     volume = volume - 5.92*mol.GetNumBonds() -14.7 * r_a - 3.8 * r_na 
#     target_mol.SetProp('VOLUME_GUEST', '{0:.4f}'.format(volume))
#     return volume

def merge_CB_guest(mol):
    """Merges the CB7 structure taken from a saved file with the structure of the guest as computed
    The argments are a mol object for the guest structure and two flags gen_ETKDG and gen_no_chiral that indicate how and if the guest structure was converged"""
    # CB_HOST = Chem.SDMolSupplier('../CB_candidate.sdf', removeHs=False)[0]
    CB_HOST = Chem.MolFromMolBlock(get_CB_BLOCK(), removeHs=False)
    complex_cb_guest = Chem.CombineMols(CB_HOST, mol)
    Chem.GetSSSR(complex_cb_guest) # Dirty fix to avoid an error
    # ff = AllChem.MMFFGetMoleculeForceField(complex_cb_guest, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(complex_cb_guest), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
    # converged = ff.Minimize(n_steps, tol, tol)
    # E_complex = ff.CalcEnergy()
    # converged = converge(complex_cb_guest)
    return complex_cb_guest

# def add_noise(matrix):
#     f
def align_xyz(list_atoms, coord_matrix):
    """This method uses a PCA on the atoms of the molecules as represented by the nuclear coordinates as points 
    and their electron cloud as a sphere of points around those nuclear coordinates
    The function returnCircleAsTuple is imported via a pybind11 generated module from c++ code"""
    angle_nbs = 6 # for angle_nbs = 3 there will be points spaced by pi/3 so three sections in inclinations and six in azimutal 
    sphere_points = get_sphere_points(list_atoms, coord_matrix)
    total_matrix = np.concatenate((coord_matrix, sphere_points),axis = 0)
    pca = PCA(n_components = 3)

    transform = pca.fit_transform(total_matrix)
    transform_coord = pca.transform(coord_matrix) 

    point_cloud = zip(transform.T[1][:], transform.T[2][:])
    height = np.max(transform.T[0][:]) - np.min(transform.T[0][:])
    # rad = make_circle(point_cloud)
    rad = returnCircleAsTuple(point_cloud)

    transform_coord_centered = transform_coord.copy()
    transform_coord_centered[:,0] = transform_coord[:,0] - np.mean(transform_coord[:,0])
    transform_coord_centered[:,1] = transform_coord[:,1] - np.mean(transform_coord[:,1])
    transform_coord_centered[:,2] = transform_coord[:,2] - np.mean(transform_coord[:,2])
    ################## CRAAAP to show cylinders
    # import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    # import mpl_toolkits.mplot3d.art3d as art3d
    # top =  np.max(transform.T[0][:])
    # bottom = np.min(transform.T[0][:])
    # index_of_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06} 
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    # ax.scatter(transform_coord.T[0][:],transform_coord.T[1][:], transform_coord.T[2][:], s=[(index_of_radii[atom]*5)**2*np.pi for atom in list_atoms], c=[index_of_radii[atom] for atom in list_atoms], alpha=0.75)
    # ax.scatter(transform.T[0][:],transform.T[1][:], transform.T[2][:], s = 2, alpha=0.75)


    # n_circles = 6
    # zs = np.linspace(bottom, top, n_circles)
    # for i in xrange(n_circles):
    #     circle=plt.Circle((rad[0],rad[1]),rad[2], fill=(i == 0 or i == n_circles-1), alpha = 0.2)
    #     ax.add_patch(circle)
    #     art3d.pathpatch_2d_to_3d(circle, z=zs[i], zdir="x")
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')
    # ax.set_title('{} ({})'.format('sample', rad[2]))
    # plt.axis('equal')
    # plt.show()

    ################# CRAAAAP
    return (rad[2], height), transform_coord_centered


def get_sphere_points(list_atoms, coord_matrix):
    """Find the thinnest cylinder dimesions where the molecule could fit add clouds of points around the atoms before the PCA (pi/3 in every direction with atomic radius)
    x = r*sin(theta)*cos(theta)
    y = r*sin(theta)*sin(theta)        theta is inclination from top (0 to pi) and phi is azimuth (0 to 2pi)
    z = r*cos(theta)
    This method uses list comprehension and return all the points representing the spheres around the atoms
    """
    index_of_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06} 
    angle_nbs = 6 # for angle_nbs = 3 there will be points spaced by pi/3 so three sections in inclinations and six in azimutal 
    sphere_points = []
    for i in xrange(len(list_atoms)):
        radius = index_of_radii[list_atoms[i]]
        top_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(0)]
        bottom_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(np.pi)]
        sphere_points.append(top_sphere_point)
        sphere_points.append(bottom_sphere_point)


    new_sphere_points = [[coord_matrix[i][0]+index_of_radii[list_atoms[i]]*np.sin(inclination_angle)*np.cos(azymuth_angle), coord_matrix[i][1]+index_of_radii[list_atoms[i]]*np.sin(inclination_angle)*np.sin(azymuth_angle), coord_matrix[i][2]+index_of_radii[list_atoms[i]]*np.cos(inclination_angle)] 
    for azymuth_angle in np.linspace(0, np.pi *2, angle_nbs*2) for inclination_angle in np.linspace(0, np.pi, angle_nbs + 1)[1:-1] for i in range(len(list_atoms))]
    sphere_points.extend(new_sphere_points)
    return sphere_points

def embed_mult_confs(mol, N_CONF, i):
    """Uses basically the ETKDG embedding with random coordinates and a rms pruning distance of 1A, clearConfs is False so that the method can be called by small batches. This ensures that the method will
    return within a reasonable amount of time and throw a timeout exception if needed. Embedding is done using two threads here."""
    #return AllChem.EmbedMultipleConfs(mol,numConfs=N_CONF, useRandomCoords=True, pruneRmsThresh=1,
    #      useExpTorsionAnglePrefs=True, useBasicKnowledge=True, enforceChirality=True, ignoreSmoothingFailures=False, randomSeed=i, boxSizeMult=2.5, maxAttempts=0, clearConfs=False, numThreads=n_threads)
    # return AllChem.EmbedMultipleConfs(mol, numConfs=N_CONF, useRandomCoords=True, pruneRmsThresh=1)
    return AllChem.EmbedMultipleConfs(mol,numConfs=N_CONF, params= AllChem.ETKDG())
def get_E(mol,ids):
    """Returns the MMFF94 estimated energy for that molecule"""
    ff = AllChem.MMFFGetMoleculeForceField(mol, confId=ids, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(mol), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
    ff.Initialize()
    return ff.CalcEnergy()
def converge(mol):
    """Converges all of the molecule built-in conformations using n_threads"""
    # AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=100000, ignoreInterfragInteractions=False, nonBondedThresh=100.0, mmffVariant='MMFF94', numThreads=n_threads)
    converged = []
    for ids in range(mol.GetNumConformers()):
        AllChem.MMFFSanitizeMolecule(mol)  
        ff = AllChem.MMFFGetMoleculeForceField(mol, confId=ids, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94'), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
        ff.Initialize()
        cf = ff.Minimize(n_steps, tol, tol)
        converged.append((cf, ff.CalcEnergy()))
    return converged


def set_energy_contributions(rdkitmol, label='Guest:'):
    """
    PRE: Takes in a rdkit mol
    POST: Assigns to it energy properties of MMFF94
          - BondStretch, AngleBend, StretchBend, OopBend, Torsion, VdW, Electrostatic
    """
    mp = AllChem.MMFFGetMoleculeProperties(rdkitmol)
    for i in range(7):
        termList = [['BondStretch', False], ['AngleBend', False],
        ['StretchBend', False], ['OopBend', False], ['Torsion', False],
        ['VdW', False], ['Electrostatic', False]]
        termList[i][1] = True
        mp.SetMMFFBondTerm(termList[0][1])
        mp.SetMMFFAngleTerm(termList[1][1])
        mp.SetMMFFStretchBendTerm(termList[2][1])
        mp.SetMMFFOopTerm(termList[3][1])
        mp.SetMMFFTorsionTerm(termList[4][1])
        mp.SetMMFFVdWTerm(termList[5][1])
        mp.SetMMFFEleTerm(termList[6][1])
        ff = AllChem.MMFFGetMoleculeForceField(rdkitmol, mp)
        rdkitmol.SetProp(label+termList[i][0], '{0:12.4f}'.format(ff.CalcEnergy()))
        # print '{0:>16s} energy: {1:12.4f} kcal/mol'.format(termList[i][0],ff.CalcEnergy())
    ff = AllChem.MMFFGetMoleculeForceField(rdkitmol, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(rdkitmol, mmffVariant='MMFF94', mmffVerbosity = 0), ignoreInterfragInteractions=False, nonBondedThresh=100.0) 
    rdkitmol.SetProp(label+'Total', '{0:12.4f}'.format(ff.CalcEnergy()))
    # w=Chem.SDWriter('./sample.sdf')
    # w.write(rdkitmol)
    # w.close()
    return rdkitmol

#@timeout(60*60)
def low_energy_low_radius_conformation(SMILES, N_CONF=500, i=42):
    n_batch = 10
    #### Creates the molecule
    mol = Chem.MolFromSmiles(SMILES)
    mol = Chem.AddHs(mol)
    E_ids_list = []
    Coords_ids_dic = {}
    #### Finds the conformations of the molecules by batches so that the molecule returns
    #### Using a random seed multiple of the number of conformers for some reason (getting different conformers?)
    # for k in range(N_CONF/n_batch):
    #     conf_ids = embed_mult_confs(mol, n_batch, k*n_batch)
    #     for ids in conf_ids:
    #         print ids
    conf_ids = embed_mult_confs(mol, N_CONF, i)
    #for ids in conf_ids:
        #print ids
    #### Converge all conformations at once
    converged = converge(mol)
    if sum([x[0] for x in converged])!=0: 
        raise Exception('Not converged my boi')
    #### Loops through the conformers and computes their energy and radii
    for ids in range(0, mol.GetNumConformers()):
        #### Computes the transformed coordinates, energies and radii of the conformers
        atom_list, coord = get_atoms_coords(Chem.MolToMolBlock(mol, confId=ids))
        cyl_data, transform_coord_centered = align_xyz(atom_list, coord)
        #### Stores the data in a list
        E_ids_list.append((ids, converged[ids][1], cyl_data[0], cyl_data[1]))
        Coords_ids_dic[ids] = transform_coord_centered

    # E_ids_radius = sorted(E_ids_list, key= lambda a:a[2])
    # E_ids_energy = sorted(E_ids_list, key= lambda a:a[1])

    # min_E_conf = E_ids_energy[0]
    # min_RAD_conf = None
    # for t in E_ids_radius:
    #     if t[1]<min_E_conf[1]+20.:
    #         min_RAD_conf = t
    #         break

    return E_ids_list, Coords_ids_dic, sum([x[0] for x in converged])

# def converge_conformers(pubchem_id):
#     print 'Converge conformers for {}'.format(pubchem_id)
#     supp = Chem.SDMolSupplier('{}_confs.sdf'.format(pubchem_id), removeHs=False)    
#     w = Chem.SDWriter('{}_min_confs.sdf'.format(pubchem_id))
#     tol = 1E-12
#     N_STEPS = 20000
#     converged = []
#     for mol in supp:
#         ff = AllChem.MMFFGetMoleculeForceField(mol, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(mol), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
#         converged.append(ff.Minimize(N_STEPS, tol, tol))
#         mol.SetProp('Energy', str(ff.CalcEnergy()))
#         w.write(mol)
#     w.close()
#     return converged


# def is_mol_converged(sdfile):
#     call('{0}obminimize -n 10 -sd {1} 2>&1 | grep CONVERGED > {1}_check_converged'.format(path_232, sdfile), shell=True)
#     with open('{}_check_converged'.format(sdfile), 'rb') as r:
#         lines = r.readlines()
#         return len(lines) == 1


# def add_hs(pubchem_id):
#     """Adds the hydrogens for the mols from the pubchem sdf files"""
#     print 'Adding Hs for {}'.format(pubchem_id)
#     call('obabel {0}.sdf -O {0}.sdf -h'.format(pubchem_id).split(' '))
#     # Popen("obabel {0}.sdf -O {0}n.sdf -h".format(pubchem_id))

# def make_sdf_from_OBABEL(SMILES, pubchem_id):
#     """Creates the 3D structure, in v.2.3.2 the hydrogens should be added anyway"""
#     print 'make obabel sdf for {}'.format(pubchem_id)
#     call('{2}obabel -:{0} -O {1}.sdf --gen3d -h'.format(SMILES, pubchem_id, path_232).split(' '))

# def make_conformers_from_sdf(pubchem_id):
#     print 'Make conformers for {}'.format(pubchem_id)
#     call('obabel {0}.sdf -O {0}_confs.sdf --conformer --nconf 5 --score rmsd --writeconformers'.format(pubchem_id).split(' '))

# @timeout(300)
# def low_energy_low_radius_conformation_OBABEL(nbr, SMILES, pb_dic):
#     #### Gets the molecule structure from PUBCHEM DICS
#     try:
#         mol_block = pb_dic[nbr]
#         if mol_block == None:
#             raise KeyError('No such mol blocks my friend')
#         w = Chem.SDWriter('{}.sdf'.format(nbr))
#         w.write(Chem.MolFromMolBlock(mol_block, removeHs=False))
#         w.close()
#         add_hs(nbr)

#     except KeyError:
#         make_sdf_from_OBABEL(SMILES, nbr)

        
#     #### Create the conformers using OBABEL 232 (ON the cluster at least)
#     print 'Generating conformers for {}'.format(SMILES)
#     make_conformers_from_sdf(nbr)

#     #### Converges the conformers using OBABEL 232
#     print 'Converging conformers for {}'.format(SMILES)
#     converge_conformers(nbr)

#     #### Do the sorting using the rdkit
#     supp_PB = Chem.SDMolSupplier('{}_min_confs.sdf'.format(nbr), removeHs=False)

#     CONVERGED_MOLS_DIC = []
#     for mols in supp_PB:
#         CONVERGED_MOLS_DIC.append(mols)
    
#     E_ids_list = []
#     Coords_ids_dic = {}

#     #### Loops through the conformers and computes their energy and radii
#     for i, mol in enumerate(CONVERGED_MOLS_DIC):
#         #### Computes the transformed coordinates, energies and radii of the conformers
#         # print('The energy inside the supp is {}'.format(get_E(mol)))
#         atom_list, coord = get_atoms_coords(Chem.MolToMolBlock(mol))
#         cyl_data, transform_coord_centered = align_xyz(atom_list, coord)
#         #### Stores the data in a list
#         E_ids_list.append((i, get_E(mol), cyl_data[0], cyl_data[1]))
#         Coords_ids_dic[i] = transform_coord_centered

#     E_ids_radius = sorted(E_ids_list, key= lambda a:a[2])
#     E_ids_energy = sorted(E_ids_list, key= lambda a:a[1])
#     w = Chem.SDWriter('best_conf_{}.sdf'.format(nbr))
#     w.write(CONVERGED_MOLS_DIC[E_ids_energy[0][0]])

#     min_E_conf = E_ids_energy[0]
#     min_RAD_conf = None
#     for t in E_ids_radius:
#         if t[1]<min_E_conf[1]+20.:
#             min_RAD_conf = t
#             break

#     # call('rm {0}.sdf {0}_confs.sdf'.format(nbr).split(' '))
#     return min_E_conf, min_RAD_conf, Coords_ids_dic[min_E_conf[0]], Coords_ids_dic[min_RAD_conf[0]]


def set_miniball_data(mol, label='Guest:'):
    atom_list, atom_coords = get_atoms_coords(Chem.MolToMolBlock(mol))
    atom_coords_as_tuples = [tuple(x) for x in atom_coords]
    miniball_data = doit(atom_coords_as_tuples)
    miniball_data[-1] = miniball_data[-1]**.5
    mol.SetProp(label+'MINIBALL(cx-xy-cz-r)', ', '.join(['{0:.4f}'.format(x) for x in miniball_data]))
    return 

times = []

def get_line_number(fname):
    """
    PRE: Gives a valid file named fname
    POST: Returns the number of lines as an int
    """
    with open(fname,'rb') as r:
        for i, lines in enumerate(r):
            pass
    return i

def one_file_treatment(fname):
    # Creates the error file for this particular input file
    err_file = fname+'_ERR'
    out_file = fname+'_OUT.sdf'
    sum_file = fname+'_SUM'
    w = Chem.SDWriter(out_file)
    with open(err_file, 'wb'): pass
    with open(sum_file, 'wb'): pass
    line_total = get_line_number(fname)

    with open(fname, 'rb') as r:
        for i, lines in enumerate(r):
            # start = time.time()
            # if i==2: break
            try: # Tries to read the line and writes the line to the error file if not succeeding
                line = lines.strip().split('\t')
                s = line[0]
                mol = Chem.MolFromSmiles(s)
                if mol.GetNumHeavyAtoms() >= 100:
                    raise Exception('TOO MANY ATOMS IN {}'.format(s))
                # s = Chem.MolToSmiles(mol, isomericSmiles=True)
                # s = strip_salts(s)
                nbr = line[1]
            except Exception, e:
                with open(err_file, 'ab') as a: a.write(lines)
                logging.error(e, exc_info=True)
                continue

            ##### BUILDS THE MOL OBJECT AND OBTAINS THE BEST ENERGY AND RADIAL CONFORMATIONS
            try:
                IDS, COORDS_IDS, converged_guest =  low_energy_low_radius_conformation(s, N_CONF=200, i=42)
                sorted_guests = sorted(IDS, key= lambda a:a[1])
                LOWEST_E = sorted_guests[0][1]
                LOWEST_E_ID = sorted_guests[0][0]
            except Exception, e: # If fails, writes the line to the err file
                with open(err_file, 'ab') as a: a.write(lines)
                logging.error(e, exc_info=True)
                continue
            print '{}/{}'.format(i,line_total), s, LOWEST_E
            # if i%20==0: print i, nbr

        #     ##### Gets a dummy mol object to serve as place holder in the actual block generation
            mol = Chem.MolFromSmiles(s)
            mol = Chem.AddHs(mol)
            conf_mol_list = []
        #     #### Gets the best energy conformation shifted by half the height of the molecule
        #     upshift_min_E_conf = min_E_conf.copy()
        #     upshift_min_E_conf[:,2] = upshift_min_E_conf[:,2]+min_E[-1]/2.0    
        #     downshift_min_E_conf = min_E_conf.copy()
        #     downshift_min_E_conf[:,2] = downshift_min_E_conf[:,2]-min_E[-1]/2.0
            for confs in sorted_guests:
                temp_mol = generate_mol_from_MDL(mol, COORDS_IDS[confs[0]])
                if temp_mol == None:
                    with open(err_file, 'ab') as a: a.write(lines)
                    print 'This molecule is None?!'
                    continue

                temp_mol.SetProp('_Name', '{}-conf-{}'.format(nbr, confs[0]))
                set_miniball_data(temp_mol)
                set_energy_contributions(temp_mol)
                temp_mol.SetProp('Radius:', '{0:12.4f}'.format(confs[2]))
                temp_mol.SetProp('Converged:Guest', '{}'.format(converged_guest==0))
                temp_mol.SetProp('SMILES', s)
                w.write(temp_mol)
                conf_mol_list.append(temp_mol)


        #     mol_best_E = generate_mol_from_MDL(mol, min_E_conf)  # Mol with the conformation corresponding to the lowest MMFF94 energy among conformers sampled
        #     mol_best_RAD = generate_mol_from_MDL(mol, min_RAD_conf) # Mol with the conformation corresponding to the lowest radius among valid (low energy) conformers sampled
        # #     mol_best_E_upshift = generate_mol_from_MDL(mol, upshift_min_E_conf, fname)  # Mol with the conformation corresponding to the lowest MMFF94 energy among conformers sampled but upshifted by h/2
        # #     mol_best_E_downshift = generate_mol_from_MDL(mol, downshift_min_E_conf, fname)  # Mol with the conformation corresponding to the lowest MMFF94 energy among conformers sampled but downshifted by h/2

        #     mol_best_E.SetProp('_Name', '{}-guest-best-E'.format(nbr))
        #     mol_best_RAD.SetProp('_Name', '{}-guest-best-RAD'.format(nbr))
        #     for best_mol, radii in zip([mol_best_E, mol_best_RAD], [min_E[2], min_RAD[2]]):
        #         set_miniball_data(best_mol)
        #         set_energy_contributions(best_mol)
        #         best_mol.SetProp('Radius:', '{0:12.4f}'.format(radii))
        #         best_mol.SetProp('Converged:Guest', '{}'.format(converged_guest==0))
        #         best_mol.SetProp('SMILES', s)
        #         w.write(best_mol)

            # set_miniball_data(mol_best_E)
            # set_miniball_data(mol_best_RAD)
            # set_energy_contributions(mol_best_E)
            # set_energy_contributions(mol_best_RAD)
            # w.write(mol_best_E)
            # w.write(mol_best_RAD)

        #     ##### GETS THE RELEVANT COMPLEXES WITH THE RELEVANT GUEST CONFORMATIONS
            complex_list = None
            try:
                for j, confs in enumerate(conf_mol_list):
                    if j == 0:
                        complex_list = merge_CB_guest(confs)
                    else:
                        ids = complex_list.AddConformer(merge_CB_guest(confs).GetConformer(), assignId=True)
                        print ids 
                # inclusion_complex_best_E  = merge_CB_guest(mol_best_E)
                # inclusion_complex_best_RAD = merge_CB_guest(mol_best_RAD)
                # exclusion_complex_up_best_E = merge_CB_guest(mol_best_E_upshift)
                # exclusion_complex_down_best_E = merge_CB_guest(mol_best_E_downshift)
            except Exception, e: # If fails, writes the line to the err file
                with open(err_file, 'ab') as a: a.write(lines)
                logging.error(e, exc_info=True)
                continue
            complex_converged = converge(complex_list)
            print complex_converged, len(complex_converged)
            ##### SETTING THE PROPERTIES
            complex_binding_energies = []
            complex_ball= []
            mol_list = []
            for i in range(complex_list.GetNumConformers()):#, exclusion_complex_up_best_E, exclusion_complex_down_best_E]): 
                print i, complex_converged[i]
                complexes = [Chem.MolFromMolBlock(Chem.MolToMolBlock(complex_list, confId=i), removeHs=False), complex_converged[i][1], complex_converged[i][0]]
                complexes[0].SetProp('_Name', '{}-inclusion-complex-attempt-{}'.format(nbr, i))
                complexes[0].SetProp('PUBCHEM_NUMBER',str(nbr))
                complexes[0].SetProp('COMPLEX_ENERGY','{0:.4f}'.format(complexes[1]))
                complexes[0].SetProp('BINDING_ENERGY','{0:.4f}'.format(complexes[1]-LOWEST_E-(-1451.415064)))
                # complexes[0].SetProp('GUEST_BEST_ENERGY_ENERGY','{0:.4f}'.format(min_E[1]))
                # complexes[0].SetProp('GUEST_BEST_RADIUS_ENERGY','{0:.4f}'.format(min_RAD[1]))
                # complexes[0].SetProp('GUEST_BEST_ENERGY_RADIUS','{0:.4f}'.format(min_E[2]))
                # complexes[0].SetProp('GUEST_BEST_RADIUS_RADIUS','{0:.4f}'.format(min_RAD[2]))
                # complexes[0].SetProp('GUEST_FORMAL_CHARGE','{}'.format(Chem.GetFormalCharge(mol_best_E)))
                complexes[0].SetProp('COMPLEX_NUMBER','{}'.format(i)) # Index in the complex list to know if it was supposed to be [0:inclusion_complex_best_E, 1:inclusion_complex_best_RAD, 2:exclusion_complex_up_best_E, 3:exclusion_complex_down_best_E]
                complexes[0].SetProp('CONVERGED','GUEST:{}, COMPLEX:{}'.format(converged_guest==0, complexes[2]==0))    
                set_miniball_data(complexes[0], label='Complex:')
                set_energy_contributions(complexes[0], label='Complex:')

                # w.write(complexes[0])
                mol_list.append(complexes[0])
                complex_binding_energies.append(complexes[1]-LOWEST_E-(-1451.415064))
                complex_converged.append(complexes[2]==0)
                complex_ball.append(complexes[0].GetProp('Complex:MINIBALL(cx-xy-cz-r)').split(',')[-1])
            # w.flush()
            complex_tuple = sorted(zip(['{}'.format(x) for x in complex_converged], ['{0:.4f}'.format(x) for x in complex_binding_energies], complex_ball, [x[0] for x in IDS]), key = lambda a: float(a[1]), reverse=False)
            ### Writes the 5 best complexes
            for k in [int(x[3]) for x in complex_tuple]:
                w.write(mol_list[k])

            with open(sum_file, 'ab') as a: #The string is long and formatted only with '\t' and 'spaces'
                fstr = '{}\t'*9 + '\n'
                a.write(fstr.format(s, nbr, LOWEST_E, LOWEST_E_ID, \
                    converged_guest==0,  complex_tuple[0][0], complex_tuple[0][1], complex_tuple[0][2], complex_tuple[0][3]))

            # stop = time.time()
            # times.append(stop-start)
        w.close()
	
	# time.sleep(30)
	
	# print fname, 'tar -cvjf {0}.tar {0}_OUT.sdf {0}_ERR {0}_SUM'.format(fname).split(' ')
	
	call('tar -cvjf {0}.tar {0}_OUT.sdf {0}_ERR {0}_SUM'.format(fname).split(' '))
	call('rm {}_OUT.sdf'.format(fname).split(' '))
        
	# print '\n'.join(['{0:.2f}'.format(x) for x in times])
        # print 'average was {}s total was {}s'.format(sum(times)/len(times), sum(times))

if __name__ == '__main__':
	# p = multiprocessing.Pool(len(flist))
	# print len(flist)
	# p.map(one_file_treatment, flist)
    start = time.time()
    one_file_treatment(flist[0])
    stop = time.time()
    print 'time is {}'.format(stop-start)
