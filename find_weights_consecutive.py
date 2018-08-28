import openbabel
import glob
# import timeit
import numpy as np
# import string
import matplotlib.pyplot as plt
from matplotlib import mlab
# file_name = sys.argv[1]

VdW_Volumes = {'C':20.58, 'N':15.60,  'S':24.43, 'Cl':22.45, 'H':7.24} #VdW volume of atoms A^3, Yuan Zhao 2003
#index_of_weight = {'H' :0, 'C':1, 'N':2, 'O':3, 'F':4, 'Cl':5, 'Br':6, 'I':7, 'P':8, 'S':9, 'As':10, 'B':11, 'Si'12, 'Se':13, 'Te':14 } 
index_of_weight = {'H' :7.24 - 5.92, 'C':20.58, 'N':15.60, 'O':14.71, 'F':13.31, 'Cl':22.45, 'Br':26.52, 'I':32.52, 'P':24.43, 'S':24.43, 'As':26.52, 'B':40.48, 'Si':38.79, 'Se':28.73, 'Te':36.62} 


# Deep Inside Cucurbiturils: Physical Properties and Volumes 
CB5_vol = 68
CB6_vol = 142
CB7_vol = 242 #55% optimal fitting
CB8_vol = 367
CB9_vol = 515 # Not partical
CB10_vol = 691 # Not practical
vols = [CB5_vol, CB6_vol, CB7_vol, CB8_vol, CB9_vol, CB10_vol]
list_of_files = glob.glob('./below_13_Hvy_Atoms_volume_dump.txt_smiles_only.can') 
# print list_of_files

def dump_volumes_new(fname, prefix):
	"""Probably better to run it on files with less than 100.000 entries for some reason, on my device it always crashes around 150.000 entries for big files"""
	"""This step erases the invalid string file"""
	with open(prefix+fname+'_INVALID_STRING', 'wb'): pass
	with open(prefix+fname+'_dumps', 'wb') as w: pass
	with open(prefix+fname+'_BUGGY_STRING', 'wb'): pass
		
	with open(prefix+fname, 'rb') as f:
		for i, lines in enumerate(f):
			if i%1==0: print '{}th line processed'.format(i), lines
			if i not in [390734]:
				smi = lines.split('\t')[0]
				volume = compute_volume_from_SMI(smi)
				if volume != 'INVALID_STRING':
					with open(prefix+fname.split('.')[0]+'_dumps.can', 'ab') as w: 
						w.write('{}\t{}\n'.format(lines.strip(), volume))
				else:
					with open(prefix+fname+'_INVALID_STRING', 'ab') as ab: 
						ab.write(lines)
			else: 
				print lines + '\n' + '#'*(10)
				with open(prefix+fname+'_BUGGY_STRING', 'ab') as ab:
					ab.write(lines)


"""HUGE AMOUNT OF CRAP"""
# def dump_volumes(file_name):
# 	""" Takes a file name in input that is expected to contain one valid SMILES fomula by line
# 		The method computes the volumes and retrieves the various indexes of molecular properties by calling the add_line_to_matrix method
# 		The method then saves the result in a text file along with the original SMILES formula
# 	"""
# 	with open(file_name, 'r') as f:
# 		global w, num, below_13_Hvy_Atoms, weight_matrix, invalid_SMILES # Global variables to be modified
# 		nb_lines = 0
# 		n = 0
# 		formula_index = []

# 		for lines in f: # Counts the lines to preallocate the matrix and build a callable list formula_index
# 			nb_lines = nb_lines+1
# 			formula_index.append(str.strip(lines))
# 		print nb_lines

# 	energy_index = []
# 	molwt_index= []
# 	dim_index= []
# 	arom_index= []
# 	narom_index= []
# 	nheavy_index= []
# 	n_rotors= []
# 	n_residues= []
# 	invalid_SMILES = []
# 	 # weights of the different matrix components in the order [Carbon atomic volume, Nitrogen atomic volume, Oxygen atomic volume, ...
# 	 # Sulfur atomic volume, Chlorine atomic volume, Hydrogen atomic volume, Bond correction, Aromatic cycle correction, Non aromatic cycle correction] 
# 	# weight_matrix = np.array([20.58, 15.60, 14.71, 24.43, 22.45, 7.24, -5.92, -14.7, -3.8])

# 	with open('dumps_%s.txt'%file_name[2:], 'w') as w:
# 		with open('volumes_in_%s_.txt'%file_name[2:], 'w') as num:
# 			with open('below_13_Hvy_Atoms_volume_dump.txt','a') as below_13_Hvy_Atoms:

# 				#w.write('%s\t%s\t%s\t%s\n'%('SMILES', 'Volume', 'Weight', 'NumHvyAtoms'))
# 				map(add_line_to_matrix, formula_index) # Effectively fills the weight matrix and the various indexes

# 				print "the corrupted smiles are %s, they won't be processed any further"%str(invalid_SMILES)
# 				"""Shrink the matrix here otherwise the reshape will fail"""

# def add_line_to_matrix(line):
# 	""" Takes a string in input and extracts molecular properties to fill in various indexes and a [atoms types and number, bonds and cycles] matrix to compute the volume
# 		The variables modified are set global to be accessed from other parts of the code, the method doesn't return anything.
# 	"""
# 	#if divmod(n, nb_lines/1000+1)[1] == 0:
# 		#print "#"*10 + " %.2f pc  complete on %i lines"%(n/float(nb_lines)*100, nb_lines)
# 	global w, num, below_13_Hvy_Atoms, weight_matrix, invalid_SMILES # Global variables to be modified

# 	line = line.split('\t')[0]
# 	obConversion.ReadString(mol, line)  # Creation of the Molecule using the openbabel c++ classes translated to python
# 	formula = mol.GetSpacedFormula().split(' ') # Extract the raw formula to obtain the contribution of each atom species
# 	volume = 0
# 	while (formula[-1] =='-') or (formula[-1] == '+'):
# 		formula.pop()
# 	atoms_are_valid = [(var in index_of_weight.keys()) for var in formula[::2]]
# 	if np.all(atoms_are_valid): # Check if the SMILES is valid or not, if not we suppose the nbr of heavy atoms will be 0, which makes no sense
# 		for a in zip(formula[::2], formula[1:][::2]): # Filling the matrix part used to compute the atom species, the numbers correspond to the number of atoms of each species occuring in the molecule
# 			volume = index_of_weight[a[0]]*float(a[1]) + volume
# 		r_a = 0	
# 		r_na = 0
# 		for r in mol.GetSSSR(): # Here the aromatic and non aromatic cycles are counted. GetSSSR() returns a list of OBRings objects that possess an IsAromatic property
# 			if r.IsAromatic():
# 				r_a = r_a+1
# 			else:
# 				r_na = r_na+1
# 		volume = volume - 5.92*mol.NumBonds() -14.7 * r_a - 3.8 * r_na
# 		output_line = '%s\t%.2f\t%.2f\t%i\n'%(line, volume, mol.GetMolWt(), mol.NumHvyAtoms())
# 		w.write(output_line)
# 		if mol.NumHvyAtoms()<=13:
# 			below_13_Hvy_Atoms.write(output_line)
# 		num.write('%.2f\n'%volume)
# 	else :
# 		invalid_SMILES.append(line)
	
def draw_cumul(file_name, isVolume):
	"""Takes a file formatted as in dump_volumes and draw its volume cumulative distribution"""
	if isVolume:
		pcs_CB7 = [24.2, 48.4, 72.6, 96.8, 121.0, 145.2, 169.4, 193.6, 217.8]
		corr = 242.
		xlab = 'Volume (A^3)'
		loc = 2
	else:
		pcs_CB7 = np.linspace(2.3, 3, 5)
		corr = 2.7
		xlab = 'radius (A)'
		loc = 1
	import colorsys
	N = len(pcs_CB7)
	HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
	RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
	x = open(file_name, 'r').readlines()
	x = np.array([float(str.strip(str.split(i, '\t')[loc])) for i in x])
	print x
	fig, ax1 = plt.subplots(figsize=(8, 4))
	n_bins1 = 50
	n, bins, patches = ax1.hist(x, n_bins1, normed=True, histtype='step', cumulative=True, label = 'Cumulative' , color = 'g')
	n_bins2 = 50
	ax2 = ax1.twinx()
	n, bins, patches = ax2.hist(x, bins=n_bins2, normed=False,  histtype='step', cumulative=False, label = 'Not Cumulative', color = 'b')
	for lns in zip(pcs_CB7, RGB_tuples):
		plt.axvline(x = lns[0], label = '{}% CB7 ({})'.format(int(lns[0]/corr*100.0), lns[0]), c = lns[1])
	#ax.plot(, , 'k--', linewidth=1.5, label='Theoretical')
	plt.legend()
	ax1.grid(True)
	ax1.set_xlim(np.min(x)*0.8, np.max(x)*1.1)
	ax1.set_title('%s'%file_name.split('/')[-1])
	ax1.set_xlabel(xlab)
	ax1.set_ylabel('Percentage below max of %.2f '%np.max(x), color = 'g')
	ax2.set_ylabel('Frequency distribution of volumes (median is %.2f)'%np.median(x), color = 'b')
	
	#fig.savefig('img_%i_%i_%s.pdf'%(n_bins1, n_bins2, file_name), format = 'pdf')
	fig.savefig("%s_img_%i_%i.pdf"%(file_name, n_bins1, n_bins2), format = 'pdf')


def get_volumes_only(fname):
	"""Takes a file formatted as in dump volumes and return a file with the volumes only"""
	with open(fname, 'rb') as f:
		with open(fname+'_vols_only.txt', 'wb') as a:
			for lines in f:
					vol = lines.split('\t')[1]
					print vol
					print lines
					a.write(vol + '\n')

def compute_volume_from_SMI(SMILES):
	"""
	Returns the volume of a molecule described by a SMILES as defined in the article
	Zhao, Yuan H., Michael H. Abraham, and Andreas M. Zissimos. "Fast calculation of van der Waals volume as a sum of atomic and bond contributions and its application to drug compounds." The Journal of organic chemistry 68.19 (2003): 7368-7373.
	REQUIRES: 
	- OPENBABEL LIBRARIES
	- Python Bindings
	- string module
	"""
	volume = 0
	index_of_vols = {'H' :7.24, 'C':20.58, 'N':15.60, 'O':14.71, 'F':13.31, 'Cl':22.45, 'Br':26.52, 'I':32.52, 'P':24.43, 'S':24.43, 'As':26.52, 'B':40.48, 'Si':38.79, 'Se':28.73, 'Te':36.62} 
	obConversion = openbabel.OBConversion()
	obConversion.SetInFormat('can')
	mol = openbabel.OBMol()
	obConversion.ReadString(mol, SMILES)
	mol.AddHydrogens() # Important otherwise OB doesn't count the covalent hydrogen bonds
	formula = mol.GetSpacedFormula().split()
	while (formula[-1] =='-') or (formula[-1] == '+'): # Strips the charges from the end of the formula
		formula.pop()
	formula_atoms = formula[::2]
	formula_atoms_nbr = formula[1:][::2]
	for at in zip(formula_atoms, formula_atoms_nbr):
		if at[0] in index_of_vols:
			volume = volume + index_of_vols[at[0]]*float(at[1])
		else: return 'INVALID_STRING'
	r_a, r_na = 0, 0
	for r in mol.GetSSSR(): # Here the aromatic and non aromatic cycles are counted. GetSSSR() returns a list of OBRings objects that possess an IsAromatic property
		if r.IsAromatic():
			r_a = r_a+1
		else:
			r_na = r_na+1
	volume = volume - 5.92*mol.NumBonds() -14.7 * r_a - 3.8 * r_na 
	return volume

if __name__ == "__main__":
	

	# prefix = '/home/hugues/Documents/mail/cb-hosts/'
	print compute_volume_from_SMI('c1ccccc1CC2CC2')
	# flist = glob.glob(prefix+'NO_*')
	# flist = [elems.split("/")[-1] for elems in flist]
	# for fnames in flist:
		# dump_volumes_new(fnames, prefix)
		# print fnames
	# prefix = '/home/hugues/Documents/mail/filter/volumes/proper_pubchem_desirable_volumes/'
	# flist = glob.glob(prefix + 'NO_UNDESIRABLE_ONLY_HETHET_ONLY_DESIRABLE_245628.can_dumps')
	# for fs in flist:
	# 	draw_cumul(fs, True)
	# 	print fs
	# print flist
