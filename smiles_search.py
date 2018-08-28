# import openbabel
# import pybel
import os
import pickle
# from file_operations import make_sample_file, file_len

n_sample = 300
# file_name = 'below_13_Hvy_Atoms_volume_dump.txt_mols_vol_108.90_to_145.20_A.txt'
# list_to_shake = ""
# file_smi = 'below_13_Hvy_Atoms_volume_dump.txt_mols_vol_108.90_to_145.20_A.txtsmiles_only.smi'
# f = open(file_name, 'r')
################### Command Line Version
# babel -ismi below_13_Hvy_Atoms_volume_dump.txt_mols_vol_108.90_to_145.20_A.txtsmiles_only.smi -O ONLY_C_Halogens.can --filter "s='[C,c][Br,F,I,Cl]' s!='[N,n,O,o,S,s,Te,Si,B,P,p,As,Se]'"
# babel -ismi below_13_Hvy_Atoms_volume_dump.txt_mols_vol_108.90_to_145.20_A.txtsmiles_only.smi -O ONLY_C_WITH_NOT_ONLY__HALOGENS.can --filter "s='[C,c][N,n,O,o,S,s,Te,Si,B,P,p,As,Se]'"
# babel -ismi below_13_Hvy_Atoms_volume_dump.txt_mols_vol_108.90_to_145.20_A.txtsmiles_only.smi -O ONLY_C_WITH_NOSPTeSe.can --filter "s='[C,c][N,n,O,o,S,s,Te,P,p,Se]'"

class Compound():
	"""Defines a compound object with attributes such as a smi_id, its cylinder data and its openbabel optimized geometry or its approximate volume as computed with the linear formula"""
	def __init__(self, smi_id = '', cylinder_height = 0, cylinder_circle = tuple, coord_xyz = '', approx_volume = 0, id_number_pubchem = 0, nb_heavy_atoms = 0):
		self.smi_id = smi_id 
		self.cylinder_height = cylinder_height
		self.cylinder_circle = cylinder_circle
		self.coord_xyz = coord_xyz
		self.approx_volume = approx_volume
		self.id_number_pubchem = id_number_pubchem
		self.nb_heavy_atoms = nb_heavy_atoms
	def __repr__(self):
		return "smi_id = %s, cylinder_height = %s, cylinder_circle = %s, coord_xyz = %s, approx_volume = %s, id_number_pubchem = %s, nb_heavy_atoms = %s"%(self.smi_id, self.cylinder_height, str(self.cylinder_circle), self.coord_xyz, self.approx_volume, self.id_number_pubchem, self.nb_heavy_atoms)

def make_fast_index(file_string = str):
	"""MOSTLY USELESS, NOT WORKING"""
	conv=openbabel.OBConversion()
	conv.OpenInAndOutFiles(file_smi,"index.fs")
	conv.SetInAndOutFormats("smi","fs")
	conv.Convert()
	conv.CloseOutFile()

def retrieve_from_index(input_name, output_name, filter_string):
	"""Applies the filter contained in filter_string to input_name file and creates (or overwrites the file output_name)"""
	os.system("obabel %s -O %s --filter %s"%(input_name, output_name, filter_string))

def convert_to_svg(input_name = str):
	"""Applies the filter contained in filter_string to input_name file and creates (or overwrites the file output_name)"""
	#print input_name[:-3]
	os.system("obabel %s -O %s"%(input_name, input_name[:-4]+'_{0}_of_{1}.svg'.format(n_sample, file_len('_'.join(input_name.split('_')[:-1])))))

def convert_to_xyz(input_name = str):
	"""Applies the filter contained in filter_string to input_name file and creates (or overwrites the file output_name)"""
	#print input_name[:-3]
	os.system("obabel %s -O %s --gen3D"%(input_name, input_name[:-4]+'.xyz'))

def make_unique(input_name, output_name):
	"""removes all duplicates and stereochemical information from the file"""
	os.system("obabel {} -O {} --unique /nostereo".format(input_name, output_name))

def convert_smi_to_can(input_name):
	"""Converts the entries of input_name in SMI to another file in CAN format"""
	os.system("obabel %s -O %s"%(input_name, input_name[:-4]+'.can'))

def make_filter_string(filter_list = list, include = True):
	"""Makes a filter string from a list of SMARTS PATTERNS"""
	if include:
		return "\"s='[\$(" + "),\$(".join(filter_list) +")]'\""
	else:
		return "\"s!='" + "' s!='".join(filter_list) +"'\""

def swap_num_smi(fname):
	with open(fname, 'rb') as f:
		with open(fname + 'swap', 'wb') as wb:
			for lines in f:
				line = lines.strip().split('\t')
				wb.write('{}\t{}\n'.format(line[1], line[0]))


def match_can_file_with_dic(dictionary_name, smi_file, outfile):
	"""Matches the SMILES entries from the smi_file with the dictionary entries and writes it to the out_file"""
	dic = {}
	with open(outfile, 'wb'): pass
	with open(dictionary_name, 'rb') as r:
		dic = pickle.load(r)
	with open(smi_file, 'r') as f:
		for lines in f: 
			line = lines.strip().split('\t')[0]
			#print 'Treating %ith compound in the smi file'%i
			if line not in dic:
				with open(outfile, 'ab') as w:
					w.write(lines)
					

def get_CB7():
	"""Returns the list of atoms and the positions as an array of the CB7 molecule, the inner radius is like 2.7 Angstroms"""

	cb7_str = """
	N        1.7467290000      3.0510900000     -4.3552700000                 
	C        1.5386000000      4.4076100000     -3.8689600000                 
	C       -0.0264900000      4.5846900000     -3.9387000000                 
	N       -0.4879800000      3.3056200000     -4.4581200000                 
	C        0.5586600000      2.4250600000     -4.7377500000                 
	N        1.8449500000      4.6676200000     -2.4699000000                 
	C        0.7155690000      5.0017800000     -1.7198200000                 
	N       -0.3906700000      4.9169900000     -2.5683200000                 
	O        0.4611800000      1.3332780000     -5.3082500000                 
	O        0.7080100000      5.3975100000     -0.5491480000                 
	N       -2.8196420000      2.7569500000     -3.8751600000                 
	C       -3.4866710000      3.7687700000     -3.0694800000                 
	C       -4.6145300000      2.9578700000     -2.3244900000                 
	N       -4.4252710000      1.5964200000     -2.8054000000                 
	C       -3.3986740000      1.4928300000     -3.7458700000                 
	N       -2.7224090000      4.3630600000     -1.9815400000                 
	C       -3.2502190000      4.0641900000     -0.7236600000                 
	N       -4.3393110000      3.2119000000     -0.9179360000                 
	O       -3.1191610000      0.4917170000     -4.4142100000                 
	O       -2.8856230000      4.5469100000      0.3540680000                 
	N        5.2614600000     -1.3809700000     -0.7897190000                 
	C        5.8213900000     -1.0076000000      0.5013800000                 
	C        5.2704690000     -2.1196300000      1.4730400000                 
	N        4.4724890000     -2.9680600000      0.5993700000                 
	C        4.5038990000     -2.5524900000     -0.7337860000                 
	N        5.3454990000      0.2308570000      1.1011000000                 
	C        4.6427100000      0.0216610000      2.2895800000                 
	N        4.5601500000     -1.3576800000      2.4902300000                 
	O        4.0293690000     -3.1613000000     -1.6986500000                 
	O        4.2523090000      0.8973510000      3.0689400000                 
	N       -1.8721010000     -4.4615300000      2.1911400000                 
	C       -2.4715060000     -3.9456000000      3.4140800000                 
	C       -3.8023090000     -3.2663300000      2.9107700000                 
	N       -3.7782690000     -3.5018600000      1.4744200000                 
	C       -2.6633280000     -4.2348500000      1.0631600000                 
	N       -1.7752000000     -2.8615300000      4.0916300000                 
	C       -2.4994710000     -1.6676200000      4.0912400000                 
	N       -3.6697380000     -1.8880600000      3.3624500000                 
	O       -2.4640220000     -4.6852600000     -0.0701770000                 
	O       -2.2080870000     -0.6373840000      4.7086300000                 
	N       -5.2811500000     -0.1475160000     -1.2878070000                 
	C       -5.8723310000      0.3687330000     -0.0615520000                 
	C       -5.7211200000     -0.8315900000      0.9493000000                 
	N       -5.0623200000     -1.8599900000      0.1571660000                 
	C       -4.8431100000     -1.4675500000     -1.1647810000                 
	N       -5.1780810000      1.4615900000      0.6038500000                 
	C       -4.6824810000      1.1030500000      1.8595100000                 
	N       -4.9644100000     -0.2512650000      2.0495100000                 
	O       -4.4304000000     -2.1878400000     -2.0801100000                 
	O       -4.1776900000      1.8662900000      2.6900000000                 
	C        3.9259100000     -4.2591900000      0.9779700000                 
	C       -1.7022400000      5.3831100000     -2.1534500000                 
	C        4.0595490000     -1.9136400000      3.7342300000                 
	C       -1.8418300000      3.0395400000     -4.9099100000                 
	N        3.8324190000      1.9909700000     -3.5761600000                 
	C        4.6289000000      2.8291600000     -2.6919000000                 
	C        5.4585200000      1.7861700000     -1.8497900000                 
	N        5.0049990000      0.5009030000     -2.3620300000                 
	C        4.0845300000      0.6285970000     -3.4042600000                 
	N        3.9182590000      3.5944300000     -1.6774000000                 
	C        4.2399390000      3.1984900000     -0.3770780000                 
	N        5.1127690000      2.1129700000     -0.4743380000                 
	O        3.6522200000     -0.2883010000     -4.1108100000                 
	O        3.8924200000      3.7648800000      0.6646400000                 
	C        3.0414200000      2.4836700000     -4.6892900000                 
	C        3.1860300000      4.8222500000     -1.9346500000                 
	C        5.7704800000      1.5594800000      0.6955400000                 
	C        5.6302100000     -0.7736800000     -2.0559800000                 
	C       -5.3593000000      0.5071780000     -2.5810000000                 
	C       -5.2230410000      2.8486600000      0.1754520000                 
	C       -4.8885000000     -3.2397200000      0.5761900000                 
	C       -4.7366500000     -0.9032000000      3.3269500000                 
	N        0.6725700000     -3.1373700000      4.1984900000                 
	C        1.1648900000     -4.3578700000      3.5766900000                 
	C        2.6539090000     -4.0027400000      3.2014900000                 
	N        2.7975600000     -2.6274600000      3.6567300000                 
	C        1.6461290000     -2.1393100000      4.2781100000                 
	N        0.5763100000     -4.7371600000      2.2995900000                 
	C        1.4947900000     -4.7013800000      1.2481000000                 
	N        2.7055500000     -4.2348400000      1.7651500000                 
	O        1.5408400000     -1.0624500000      4.8752500000                 
	O        1.3003500000     -5.1039300000      0.0958950000                 
	C       -0.5882600000     -3.0278700000      4.9109900000                 
	C       -0.7285100000     -5.3563600000      2.1473600000                 
	H        2.0833600000      5.1155500000     -4.5048100000                 
	H       -0.3415400000      5.3921000000     -4.6104700000                 
	H       -3.8787030000      4.5586100000     -3.7213700000                 
	H       -5.6275600000      3.2965900000     -2.5722200000                 
	H        6.9158700000     -0.9908000000      0.4384120000                 
	H        6.0621490000     -2.7142800000      1.9444000000                 
	H       -2.6450060000     -4.7690100000      4.1171800000                 
	H       -4.7066310000     -3.7110500000      3.3430300000                 
	H       -6.9158700000      0.6515360000     -0.2450230000                 
	H       -6.6814700000     -1.2094400000      1.3199700000                 
	H        4.6916400000     -4.8176600000      1.5314500000                 
	H        3.6896190000     -4.7811200000      0.0502540000                 
	H       -1.5716500000      5.8609200000     -1.1821140000                 
	H       -2.0623090000      6.1195400000     -2.8825300000                 
	H        3.9017700000     -1.0747100000      4.4129900000                 
	H        4.8201300000     -2.5869700000      4.1502400000                 
	H       -1.7916300000      2.1521600000     -5.5417200000                 
	H       -2.1899240000      3.8957700000     -5.5017600000                 
	H        5.2573390000      3.5014100000     -3.2888400000                 
	H        6.5418900000      1.8773400000     -1.9912400000                 
	H        2.8494090000      1.6300500000     -5.3402100000                 
	H        3.6263590000      3.2379000000     -5.2305600000                 
	H        3.0854500000      5.3334200000     -0.9767350000                 
	H        3.7715500000      5.4427200000     -2.6247500000                 
	H        5.5455390000      2.2278800000      1.5272300000                 
	H        6.8530600000      1.5430800000      0.5169300000                 
	H        5.3223300000     -1.4691700000     -2.8379700000                 
	H        6.7197390000     -0.6465500000     -2.0820300000                 
	H       -5.1386700000     -0.2535490000     -3.3306000000                 
	H       -6.3802200000      0.8832390000     -2.7255200000                 
	H       -4.9230000000      3.4505500000      1.0336300000                 
	H       -6.2538810000      3.0963500000     -0.1075790000                 
	H       -4.7040810000     -3.8198900000     -0.3281800000                 
	H       -5.8155000000     -3.5817200000      1.0530100000                 
	H       -5.6694300000     -1.3828900000      3.6483800000                 
	H       -4.4639200000     -0.1195970000      4.0343100000                 
	H        1.0819400000     -5.1924100000      4.2828600000                 
	H        3.3918600000     -4.6376900000      3.7059600000                 
	H       -0.5156210000     -2.1411400000      5.5417200000                 
	H       -0.7147000000     -3.9186600000      5.5402600000                 
	H       -0.8503200000     -6.1195400000      2.9258400000                 
	H       -0.7384610000     -5.8283000000      1.1644700000     
	"""

	coord_matrix = np.reshape(np.array(cb7_str.split()), (126,4))
	list_atoms = [elems[0] for elems in coord_matrix]
	coord = [elems[1:] for elems in coord_matrix]
	coord_float = np.zeros([126,3])
	for i in xrange(126):
		for j in xrange(3):
			coord_float[i][j] = float(coord[i][j])

	return list_atoms, coord_float

unstable_list_linux = {
	'[\$([CX2](=C)=C)]':'Allenes',
	'[CX3](=[OX1])[OX2][CX3](=[OX1])':'Anhydrides',
	'[CX4]([NX3])([NX3])':'Aminal_cyclic_or_not',
	#'[CX4H]([OX2H])([OX2H])',#Hemiacetal gets hit by the gem diol expression
	'[CX4]([OX2H])([OX2H])':'Gem-diols',
	'[OR0X2H0]-[CR0X4H]-[OR0X2H0]':'Acyclic_acetals', #would hit also the gem-diols without the acyclic constrain 75
	#idem as up #Acetal chains
	#Including primary amine and (ketone or aldehyde)
	"[\$([CX3]([#6])[#6]),\$([CX3H][#6])]=[\$([NX2][#6]),\$([NX2H])]":'Imines_substituted_or_not',# unstable hydrolysis 3361
	"[#7,#8]-[CR0]=[NR0]-[CR0]=[#7,#8]":'Odd_formula',
	"[OX2H][#6X3]=[#6]":'Enols',
	'[NX3][CX3]=[CX3]':'Enamines',
	'[CR0X3](=O)([OR0])[OR0]':'Acyclic_carbonates',
	'[CX3](=[OX1])(O)O':'Carbonic_acid_and_derivatives',
	'[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]':'Carbamic_acids_and_derivatives',
	'[CX3](=O)[CX4][CX3](=O)(O)':'Beta-carboxylic_acids_and_derivatives',
	#Heavily twisted amides
	'[N,n,O,o,P,p,S,s,Se,Te][F,Br,Cl,I]':'Heteroatom-halogen_bonds',
	'[#6][CX3](=O)[OX2H0][#6]':'Esthers',# also hits anhydride 3308
	'[CX3](=[OX1])[F,Cl,Br,I]':'Acid_Halides',
	#Sulfonyl halides
	#Halopyrimidines
	'[OX2R1]1[CX4R1][CX4R1]1':'Epoxies',
	'[NX3R1]1[CX4R1][CX4R1]1':'Azyridines',
	'[CX3H1](=O)[#6]':'Aldehydes',# unstable in water and in acidic conditions 4097
	'S([#6])[CX3](=O)[#6]':'Thioesters',# unstable in water 375
	'[SX4](=O)(=O)([OX2H0])[#6]':'Sulfonate_Esters',# 441
	'[PX4](=O)(-O)([OX2H0])[#6]':'Phosphonate_Esters', 
	# #'[cR1]1[cR1][cR1][cR1][cR1][cR1]1', # unfused benzene 1749
	#'[aR1]1[aR1][aR1][aR1][aR1][aR1]1', # 6 membered aromatic rings 15033
	#'[aR1]1[aR1][aR1][aR1][aR1]1' # 5 members aromatic rings 25379
	'[a]':'Any_aromatic_atom'
	}
unstable_dic_windows = { #126424 if apply all remains 42333
	'[$([CX2](=C)=C)]':'Allenes',
	'[CX3](=[OX1])[OX2][CX3](=[OX1])':'Anhydrides',
	'[CX4]([NX3])([NX3])':'Aminal_cyclic_or_not',
	#'[CX4H]([OX2H])([OX2H])',#Hemiacetal gets hit by the gem diol expression
	'[CX4]([OX2H])([OX2H])':'Gem-diols',
	'[OR0X2H0]-[CR0X4H]-[OR0X2H0]':'Acyclic_acetals', #would hit also the gem-diols without the acyclic constrain 75
	#idem as up #Acetal chains
	#Including primary amine and (ketone or aldehyde)
	"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]":'Imines_substituted_or_not',# unstable hydrolysis 3361
	"[#7,#8]-[CR0]=[NR0]-[CR0]=[#7,#8]":'Odd_formula',
	"[OX2H][#6X3]=[#6]":'Enols',
	'[NX3][CX3]=[CX3]':'Enamines',
	'[CR0X3](=O)([OR0])[OR0]':'Acyclic_carbonates',
	'[CX3](=[OX1])(O)O':'Carbonic_acid_and_derivatives',
	'[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]':'Carbamic_acids_and_derivatives',
	'[CX3](=O)[CX4][CX3](=O)(O)':'Beta-carboxylic_acids_and_derivatives',
	#Heavily twisted amides
	'[N,n,O,o,P,p,S,s,Se,Te][F,Br,Cl,I]':'Heteroatom-halogen_bonds',
	'[#6][CX3](=O)[OX2H0][#6]':'Esthers',# also hits anhydride 3308
	'[CX3](=[OX1])[F,Cl,Br,I]':'Acid_Halides',
	#Sulfonyl halides
	#Halopyrimidines
	'[OX2R1]1[CX4R1][CX4R1]1':'Epoxies',
	'[NX3R1]1[CX4R1][CX4R1]1':'Azyridines',
	'[CX3H1](=O)[#6]':'Aldehydes',# unstable in water and in acidic conditions 4097
	'S([#6])[CX3](=O)[#6]':'Thioesters',# unstable in water 375
	'[SX4](=O)(=O)([OX2H0])[#6]':'Sulfonate_Esters',# 441
	'[PX4](=O)(-O)([OX2H0])[#6]':'Phosphonate_Esters', 
	# #'[cR1]1[cR1][cR1][cR1][cR1][cR1]1', # unfused benzene 1749
	#'[aR1]1[aR1][aR1][aR1][aR1][aR1]1', # 6 membered aromatic rings 15033
	#'[aR1]1[aR1][aR1][aR1][aR1]1' # 5 members aromatic rings 25379
	'[a]':'Any_aromatic_atom'
	}

### Includes easy protonation and bonds between heavy atoms if apply from those passing the unstable test it remains 46072 => 17637
desirable_list_linux = {
	# '[N,n,O,o,P,p,S,s,Se,Te]~[N,n,O,o,P,p,S,s,Se,Te]':'Bonds_between_heteroatoms',
	'[\$([NX3](=O)=O),\$([NX3+](=O)[O-])][!#8]': 'Nitro_groups',
	'[#6][OX2H]':'Alcohols',
	'[NX1]#[CX2]':'Nitriles',
	'[#6][CX3](=O)[#6]':'Ketones',
	'[#16X2H]':'Thiols',
	"[OX2H][cX3][c]":'Phenol-like',
	'[#6][CX3](=O)[OX2H0][#6]':'Esters',
	'[#16X2H0]':'Sulfide',
	'[CX3](=O)[OX2H1]':'Carboxilic_acids',
	'[OD2]([#6])[#6]':'Ethers',
	'[NX3][CX3](=[OX1])[#6]':'Amides',
	'[NX3][cc]':'Anilines',
	# #Pyridines
	# #Imidazole
	'[NX3;H2,H1;!\$(NC=O)]':'Primary_or_secondary_Amines'
	}
desirable_dic_windows = {
	# '[N,n,O,o,P,p,S,s,Se,Te]~[N,n,O,o,P,p,S,s,Se,Te]':'Bonds_between_heteroatoms',
	'[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]': 'Nitro_groups',
	'[#6][OX2H]':'Alcohols',
	'[NX1]#[CX2]':'Nitriles',
	'[#6][CX3](=O)[#6]':'Ketones',
	'[#16X2H]':'Thiols',
	"[OX2H][cX3][c]":'Phenol-like',
	'[#6][CX3](=O)[OX2H0][#6]':'Esters',
	'[#16X2H0]':'Sulfide',
	'[CX3](=O)[OX2H1]':'Carboxilic_acids',
	'[OD2]([#6])[#6]':'Ethers',
	'[NX3][CX3](=[OX1])[#6]':'Amides',
	'[NX3][cc]':'Anilines',
	# #Pyridines
	# #Imidazole
	'[NX3;H2,H1;!$(NC=O)]':'Primary_or_secondary_Amines'
	}
looks_like_dic = {'C1CC2CCC1N=N2':'DBO1_like', '[#6]1[#6][#6]2[#6][#6][#6]1[#7]=[#7]2':'DBO2_like', '[CX4]1[CX4][CX3]2[CX4][CX4][CX4]1[NX2]=[NX2]2':'DBO3_like', 'C1CC2CC1N=N2':'DBH_like', 'C1=COC=C1':'Furan_like', 'C1CCC(=O)C1':'Cyclopentanone_like'}
looks_like_list = looks_like_dic.keys()


desirable_list_windows = desirable_dic_windows.keys()
unstable_list_windows = unstable_dic_windows.keys()

# for key in unstable_dic_windows:
# 	print key, unstable_dic_windows[key]


# convert_to_xyz(prefix+'No_AROM.can')
# non_arom_string = make_filter_string(['[a]'], False)
# retrieve_from_index(input_name = prefix+'below_13_Hvy_Atoms_volume_dump.txt_smiles_only.can', output_name=prefix+'No_AROM.can', filter_string = non_arom_string)
# print non_arom_string

print make_filter_string(unstable_list_linux.keys(), False)

def make_all_samples():
	exclude_string_linux = make_filter_string(unstable_list_linux, False)
	exclude_string_windows = make_filter_string(unstable_list_windows, False)
	include_string_windows = make_filter_string(desirable_list_windows, True)
	include_string_linux = make_filter_string(desirable_list_linux, True)
	look_alike_string = make_filter_string(looks_like_list, True)
	for items in desirable_list_windows:
		fil = make_filter_string([items], True)
		print fil
		outf = prefix+'ONLY_%s.can'%desirable_dic_windows[items]
		retrieve_from_index(input_name = prefix+'below_13_Hvy_Atoms_volume_dump.txt_smiles_only.can', output_name=outf, filter_string = fil)
		make_sample_file(outf, n_sample)
		convert_to_svg(outf+'_sample.can')

	for items in unstable_list_windows:
		fil = make_filter_string([items], True)
		print fil
		outf = prefix+'ONLY_%s.can'%unstable_dic_windows[items]
		retrieve_from_index(input_name = prefix+'below_13_Hvy_Atoms_volume_dump.txt_smiles_only.can', output_name=outf, filter_string = fil)
		make_sample_file(outf, n_sample)
		convert_to_svg(outf+'_sample.can')

	for items in looks_like_list:
		fil = make_filter_string([items], True)
		print fil
		outf = prefix+'ONLY_%s.can'%looks_like_dic[items]
		retrieve_from_index(input_name = prefix+'below_13_Hvy_Atoms_volume_dump.txt_smiles_only.can', output_name=outf, filter_string = fil)
		make_sample_file(outf, n_sample)
		convert_to_svg(outf+'_sample.can')


def make_only_final_set():
	"""	First part removes all that is not desired and stores it to outf_not_desired
		Second part is keep only those that have Bonds_between_heteroatoms
		third part is making small files for each of the desired species
		fourth part applies all the desirable filters simultaneously
	"""
	out_ns = prefix + 'merged_13_pubchem_NS.can'
	make_unique(prefix + 'merged_13_pubchem_properswap.smi', out_ns)



	outf_not_desired = prefix + 'NO_UNDESIRABLE.can'
	exclude_string = make_filter_string(unstable_list_linux.keys(), False)
	retrieve_from_index(input_name = out_ns, output_name=outf_not_desired, filter_string = exclude_string)

	"""Second part """
	outf_not_desired_only_HetHet_bond = prefix + 'NO_UNDESIRABLE_ONLY_HETHET.can'
	include_string = make_filter_string(['[N,n,O,o,P,p,S,s,Se,Te]~[N,n,O,o,P,p,S,s,Se,Te]'], True)
	print include_string
	retrieve_from_index(input_name = outf_not_desired, output_name = outf_not_desired_only_HetHet_bond, filter_string = include_string)


	"""Third part"""
	for items in desirable_list_linux.keys():
		outf = prefix + 'NO_UNDESIRABLE_ONLY_HETHET_{}.can'.format('ONLY_'+desirable_list_linux[items])
		include_string = make_filter_string([items], True)
		retrieve_from_index(input_name = outf_not_desired_only_HetHet_bond, output_name=outf, filter_string = include_string)
		make_sample_file(outf, n_sample)
		convert_to_svg(outf+'_sample.can')

	"""Fourth part"""
	retrieve_from_index(input_name = prefix + 'NO_UNDESIRABLE_ONLY_HETHET.can', output_name=prefix+"NO_UNDESIRABLE_ONLY_HETHET_ONLY_DESIRABLE.can", filter_string = make_filter_string(desirable_list_linux.keys(), True))

# retrieve_from_index(input_name = prefix + 'NO_UNDESIRABLE_ONLY_HETHET.can', output_name=prefix+"NO_UNDESIRABLE_ONLY_HETHET_ONLY_DESIRABLE.can", filter_string = make_filter_string(desirable_list_windows, True))
# prefix = './filter/Cylinders/'
# match_can_file_with_dic(prefix+'xyz_dict.txt',prefix+'NO_UNDESIRABLE_ONLY_HETHET_ONLY_DESIRABLE_245628_dumps.can', prefix+'not_present.can')

# make_only_final_set()
# swap_num_smi(prefix+'merged_13_pubchem_proper')
# make_unique('/home/hugues/Documents/mail/filter/atm_nbr_split/pubchem_total.smi_5_atoms.can')
# make_only_final_set()
# make_sample_file(prefix+'test_DBO_3Strings_final_set.can', 300)
# convert_to_svg(prefix+'test_DBO_3Strings_final_set.can_sample.can')

#convert_to_smiles_only('Compound_036200001_036225000.sdf.gz.smi')
#convert_to_svg('ONLY_NO_UNDESIRABLE_ONLY_DESIRABLE.can')
#convert_to_xyz('ONLY_CNOSPTeSe_WITH_HETEROATOM_BONDS.can')

#convert_smi_to_can('below_13_Hvy_Atoms_volume_dump.txt_smiles_only.smi')


#match_can_file_with_dic('dic_hvy_binary_dic', 'ONLY_NO_UNDESIRABLE_ONLY_DESIRABLE.can', 'ONLY_NO_UNDESIRABLE_ONLY_DESIRABLE_FITTING_CYL.can')
#convert_to_svg('ONLY_NO_UNDESIRABLE_ONLY_DESIRABLE_FITTING_CYL.can')
