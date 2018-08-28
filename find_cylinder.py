import numpy as np
import scipy
from sklearn.decomposition import PCA
from smallestenclosingcircle import make_circle
import pickle
import cPickle
import sys
import collections as c 
# import file_operations
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

# file_name = sys.argv[1]

compound = c.namedtuple('compound', ['atm_list', 'xyz_array', 'atm_nbr', 'volume', 'pubchem_nbr'])

def get_best_conformation(smiles = str):
	import openbabel
	mol = openbabel.OBMol()
	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats("smi", "xyz")
	obConversion.ReadString(mol, smiles)
	mol.AddHydrogens()
	builder = openbabel.OBBuilder()
	builder.Build(mol)
	ff = openbabel.OBForceField.FindForceField("MMFF94")
	ff.Setup(mol)
	ff.SteepestDescent(250)
	ff.WeightedRotorSearch(200,25)
	ff.ConjugateGradients(200)
	raw_xyz = obConversion.WriteString(mol)
	#print raw_xyw
	raw_string = str.split(raw_xyz,'\n')
	return raw_xyz, ff.Energy()

def compose_matrix_from_raw_xyz(xyz_raw):
	# import numpy as np
	raw_string = xyz_raw
	list_atoms = []
	coord_matrix = []
	# print raw_string
	# print 'new string'
	for lines in raw_string.split('\n')[2:-1]:
		lines = lines.split()
		# lines = " ".join(lines.split()).split(' ')
		if ('-1.#IND0' in lines or '-nan' in lines or len(lines)!=4):
			return False
		list_atoms.append(lines[0])
		coord_matrix.append(lines[1:])
	return list_atoms, np.array(coord_matrix, dtype=float)

def find_best_cylinder(list_atoms, coord_matrix):
	"""Find the thinnest cylinder dimesions where the molecule could fit add clouds of points around the atoms before the PCA (pi/3 in every direction with atomic radius)
	x = r*sin(theta)*cos(theta)
	y = r*sin(theta)*sin(theta)        theta is inclination from top (0 to pi) and phi is azimuth (0 to 2pi)
	z = r*cos(theta)
	"""
	# import numpy as np
	# import scipy
	# from sklearn.decomposition import PCA
	# from smallestenclosingcircle import make_circle
	index_of_volumes = {'H' :7.24, 'C':20.58, 'N':15.60, 'O':14.71, 'F':13.31, 'Cl':22.45, 'Br':26.52, 'I':32.52, 'P':24.43, 'S':24.43, 'As':26.52, 'B':40.48, 'Si':38.79, 'Se':28.73, 'Te':36.62} 
	index_of_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06} 
	angle_nbs = 6 # for angle_nbs = 3 there will be points spaced by pi/3 so three sections in inclinations and six in azimutal 
	sphere_points = []
	for i in xrange(len(list_atoms)):
		radius = index_of_radii[list_atoms[i]]
		top_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(0)]
		bottom_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(np.pi)]
		sphere_points.append(top_sphere_point)
		sphere_points.append(bottom_sphere_point)
		for inclination_angle in np.linspace(0, np.pi, angle_nbs + 1)[1:-1]:
			for azymuth_angle in np.linspace(0, np.pi *2, angle_nbs*2):
				new_sphere_point = [coord_matrix[i][0]+radius*np.sin(inclination_angle)*np.cos(azymuth_angle), coord_matrix[i][1]+radius*np.sin(inclination_angle)*np.sin(azymuth_angle), coord_matrix[i][2]+radius*np.cos(inclination_angle)]
				sphere_points.append(new_sphere_point)
	# print len(sphere_points)
	# print len(coord_matrix[:])
	# print np.array(sphere_points)
	total_matrix = np.concatenate((coord_matrix, sphere_points),axis = 0)
	pca = PCA(n_components = 3)

	transform = pca.fit_transform(total_matrix)
	transform_coord = pca.transform(coord_matrix) 

	point_cloud = zip(transform.T[1][:], transform.T[2][:])
	height = np.max(transform.T[0][:]) - np.min(transform.T[0][:])
	rad = make_circle(point_cloud)


	return [rad[2], height]

def show_cylinder(smiles):
	conf, energ = get_best_conformation(smiles)
	check_for_invalid_ = compose_matrix_from_raw_xyz(conf)
	print 'this is the python plot coordinate string'
	print conf
	list_atoms, coord_matrix = check_for_invalid_

	index_of_volumes = {'H' :7.24, 'C':20.58, 'N':15.60, 'O':14.71, 'F':13.31, 'Cl':22.45, 'Br':26.52, 'I':32.52, 'P':24.43, 'S':24.43, 'As':26.52, 'B':40.48, 'Si':38.79, 'Se':28.73, 'Te':36.62} 
	index_of_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06} 
	angle_nbs = 6 # for angle_nbs = 3 there will be points spaced by pi/3 so three sections in inclinations and six in azimutal 
	sphere_points = []
	for i in xrange(len(list_atoms)):
		radius = index_of_radii[list_atoms[i]]
		top_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(0)]
		bottom_sphere_point =  [coord_matrix[i][0], coord_matrix[i][1], coord_matrix[i][2]+radius*np.cos(np.pi)]
		sphere_points.append(top_sphere_point)
		sphere_points.append(bottom_sphere_point)
		for inclination_angle in np.linspace(0, np.pi, angle_nbs + 1)[1:-1]:
			for azymuth_angle in np.linspace(0, np.pi *2, angle_nbs*2):
				new_sphere_point = [coord_matrix[i][0]+radius*np.sin(inclination_angle)*np.cos(azymuth_angle), coord_matrix[i][1]+radius*np.sin(inclination_angle)*np.sin(azymuth_angle), coord_matrix[i][2]+radius*np.cos(inclination_angle)]
				sphere_points.append(new_sphere_point)
	# print len(sphere_points)
	# print len(coord_matrix[:])
	# print np.array(sphere_points)
	total_matrix = np.concatenate((coord_matrix, sphere_points),axis = 0)
	pca = PCA(n_components = 3)

	transform = pca.fit_transform(total_matrix)
	transform_coord = pca.transform(coord_matrix) 

	point_cloud = zip(transform.T[1][:], transform.T[2][:])
	height = np.max(transform.T[0][:]) - np.min(transform.T[0][:])
	rad = make_circle(point_cloud)

	top =  np.max(transform.T[0][:])
	bottom = np.min(transform.T[0][:])
	# print height
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	# print [index_of_radii[atom]**2*np.pi for atom in list_atoms]
	
	### Not transformed coordinates scatter plot
	### ax.scatter(coord_matrix.T[0][:],coord_matrix.T[1][:], coord_matrix.T[2][:], s=[(index_of_radii[atom]*5)**2*np.pi for atom in list_atoms], c=[index_of_radii[atom] for atom in list_atoms], alpha=0.75)
	### ax.scatter(total_matrix.T[0][:],total_matrix.T[1][:], total_matrix.T[2][:], s = 2, alpha=0.75)

	#### Transformed coordinates scatter plot
	ax.scatter(transform_coord.T[0][:],transform_coord.T[1][:], transform_coord.T[2][:], s=[(index_of_radii[atom]*5)**2*np.pi for atom in list_atoms], c=[index_of_radii[atom] for atom in list_atoms], alpha=0.75)
	ax.scatter(transform.T[0][:],transform.T[1][:], transform.T[2][:], s = 2, alpha=0.75)

	### Visual indications, mostly circles to outline the perimeter of the cylinder
	
	n_circles = 6
	zs = np.linspace(bottom, top, n_circles)
	for i in xrange(n_circles):
		circle=plt.Circle((rad[0],rad[1]),rad[2], fill=(i == 0 or i == n_circles-1), alpha = 0.2)
		ax.add_patch(circle)
		art3d.pathpatch_2d_to_3d(circle, z=zs[i], zdir="x")
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	ax.set_title('{} ({})'.format(smiles, rad[2]))
	plt.axis('equal')
	plt.show()


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def write_molecules_and_cylinders(input_file):
	total = file_len(input_file)
	numline = 0
	with open(input_file, 'r') as f:
		for line in f:
			numline+=1
			percentage = 100 * numline / total
			if numline%500 == 0: print line, 'percentage %.4f pc'%percentage
			conf, energy = get_best_conformation(line.strip())
			check_for_invalid_ = compose_matrix_from_raw_xyz(conf)
			if check_for_invalid_:
				list_atoms, coord_matrix = check_for_invalid_
				rad, height = find_best_cylinder(list_atoms = list_atoms, coord_matrix = coord_matrix)
				# print rad, height
				with open(input_file+'.txt', 'ab') as w:
					w.write('%s\t%s\t%s\t%s\n%s\n'%(line, rad, height, energy, conf))
				# with open(input_file+'.txt', 'ab') as w: ###########################################################REMOVE
				# 	w.write('%s\t%s\n%s\n'%(line, energy, conf))

			else: 
				with open(input_file+'_wrong.txt', 'ab') as w_wrong:
					w_wrong.write(line + '\n')

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




def make_volume_dictionary(input_files, prefix):
	"""returns a dictionary with the canon smiles as key and a tuple (Pubchem Number, approx volume) provided a text file with dump_volumes_new 'SMI\tPUBCHEMID\tVOLUME\n'"""
	keyset = set()
	pubchem_dict = {}
	with open(prefix+input_files, 'rb') as f:
		for i, lines in enumerate(f):
			lines = lines.split('\t')
			keyset.add(lines[0])
			if i%5000==0: print '{}th line added to dic'.format(i), lines
		pubchem_dict = dict.fromkeys(keyset)
		f.seek(0)
		for j, lines in enumerate(f):
			lines = lines.split('\t')
			pubchem_dict[lines[0]] = (lines[1],float(lines[2]))
			if j%5000==0: print '{}th line added to dic'.format(j)

		with open(prefix+input_files+'_number_dict.txt', 'wb') as w:
			pickle.dump(pubchem_dict, w)
		
	return pubchem_dict




def make_cylinder_dictionary(xyzdic, prefix):
	"""Returns a dictionary with can smiles as keys and named tuples as entries

	"""
	dic_cyl = dict.fromkeys(xyzdic.keys(),[])
	with open(prefix+'smiles_wrong.txt', 'ab'): pass
	for i, key in enumerate(xyzdic.keys()):
		if i%100==0: print '{}th key processed it is {}'.format(i, key)
		# print 'the key is {}'.format(key)
		# print 'the raw string is {}'.format(xyzdic[key])
		check_for_invalid_ = compose_matrix_from_raw_xyz(xyzdic[key])
		if check_for_invalid_:
			atom_list, xyz_matrix = compose_matrix_from_raw_xyz(xyzdic[key])
			dic_cyl[key] = find_best_cylinder(atom_list, xyz_matrix)
		else:
			with open(prefix+'smiles_wrong.txt', 'ab') as w_wrong:
				w_wrong.write(key + '\n')

	with open(prefix + 'cylinder_dict.txt', 'wb') as wb:
		pickle.dump(dic_cyl, wb)


def associate_err_smi_xyz(xyzfile):
	"""returns one lists with the id tuples, the xyz strings are formatted as below
	5\n
	0010100101 0101001.00\n
	C          0.94849       -0.05712       -0.03387\n
	C          2.46766       -0.02152        0.04250\n
	N          3.14712        0.71080        0.89338\n
	N          2.45397        1.50028        1.74114\n
	N          3.11101       -0.84959       -0.87546\n
	Returns False if the xyz file has an odd number of lines (which indicates it is probably corrupted) or if the number of molecules in xyz and the smifile don't match
	"""
	xyzlist = []

	with open(xyzfile, 'rb') as r:
		for i,lines in enumerate(r):
			s = ''
			if lines.strip().isdigit():
				try:
					pbId_vol = r.next().strip()
					atm_number = lines.strip()
					# print pbId_vol
					atm_list = []
					xyz_array = np.zeros([int(atm_number), 3])
					for j in range(int(lines.strip())):
						tmp = r.next()
						tmp = tmp.strip().split()
						atm_list.append(tmp[0])
						xyz_array[j, :] = np.asarray(tmp[1:])
					tmp_tuple = compound(atm_list, xyz_array, atm_number, pbId_vol.split('\t')[1], pbId_vol.split('\t')[0] )
					xyzlist.append(tmp_tuple)
				except ValueError: 
					print 'WOW {}'.format(pbId_vol), j
					with open(prefix+'xyz_problematic_string', 'ab') as ab:
						ab.write(format(pbId_vol)+'\n')
					pass
				except StopIteration:
					print 'Unfinished business with {}'.format(pbId_vol)
	return xyzlist

def merge_dicts(dic1, dic2):
	dic3 = dic1.copy()
	dic3.update(dic2)
	return dic3

def demo():
	with open('/home/hugues/Documents/mail/filter/Cylinders/XYZ_DICT_TOT.txt', 'rb') as r:
		cyl_dic = pickle.load(r) 
	print cyl_dic['NCC12CCC(CC1)N=N2']
	show_cylinder('NCC12CCC(CC1)N=N2')
	with open('/home/hugues/Documents/mail/filter/Cylinders/XYZ_DICT_TOT.txt', 'rb') as r:
		dic_xyz = pickle.load(r)
		print 'this is the dictionary string'
		print dic_xyz['C1OCB(C1)N=NB1COCC1']

# prefix = './cylinder_enclosed/cylinders_by_13/3dCoord_data/3dCoord_data/'
prefix = './filter/Cylinders/'
flist = glob.glob(prefix+'*.xyz')
flist = [el.split('/')[-1] for el in flist]

def find_intersection():
	with open('/home/hugues/Documents/mail/filter/volumes/proper_pubchem_desirable_volumes/NO_UNDESIRABLE_ONLY_HETHET_ONLY_DESIRABLE_245628.can_dumps', 'rb') as r:
		with open('/home/hugues/Documents/mail/filter/Cylinders/CYLINDER_DICT_TOT.txt', 'rb') as rb:
			dic_cyl = pickle.load(rb)
			count = 0.
			bad_key_nbr = 0
			for line in r:
				vol = float(line.split('\t')[2])
				smi = line.split('\t')[0]
				# print vol
				if (vol > 0.8*242) and (vol<1*242.):
					try:
						if (dic_cyl[smi][0][2] < 4):
							count+=1
							print smi, dic_cyl[smi][0][2], vol
					except: 
						bad_key_nbr +=1
						print '{}th bad key {}'.format(bad_key_nbr, smi)
			print count
def get_atoms_coords(RDKIT_BLOCK):
	"""Takes as input an RDKIT BLOCK and returns a list of atoms with a numpy array containing the coordinates"""
	atm_number = int(RDKIT_BLOCK[4].split()[0])
	print '#'*20
	atm_list = []
	coords_array = np.zeros([atm_number, 3], dtype=float)
	for i, line in enumerate(RDKIT_BLOCK[5:5+atm_number]):
		coords_atm = line.split()
		atm_list.append(coords_atm[3])
		coords_array[i, :] = coords_atm[:3]
	return atm_list, coords_array


def parse_RDKIT_XYZ_FILE(infile):
	"""Parse a file of concatenated RDKIT BLOCKs to find their cylinder radius and height in the PCA sense, a similar 
	file is returned where the name of the molecule is now incremented by radius and height"""
	outfile = infile+'_OUT'
	with open(outfile, 'wb'):pass
	with open(infile, 'rb') as r:
		tmp_block = []
		for i, line in enumerate(r):
			tmp_block.append(line)
			if line.strip() == "M  END": 
				atm_list, atm_coords = get_atoms_coords(tmp_block)
				cylinder = find_best_cylinder(atm_list, atm_coords)
				tmp_block[0] = tmp_block[0].strip() + '\t{0:.2f}\t{1:.2f}\t{2}\n'.format(cylinder[0], cylinder[1], '(smi, #, A^3, A, A)')
				with open(outfile, 'ab') as ab:
					ab.write(''.join(tmp_block))
				tmp_block = []

def match_parsed_RDKIT_XYZ_FILE_with_ZINC(zinc_dic_file, list_of_parsed_files):
	out_vol_cyl_id_nbs_dic = '/'.join(list_of_parsed_files[0].split('/')[:-1])+'/DATA_DIC_STEREO.txt'
	print out_vol_cyl_id_nbs_dic
	keyset = set()
	for fname in list_of_parsed_files:
		with open(fname, 'rb') as r:
			tmp_block = []
			for i, line in enumerate(r):
				if i%100000==0: print '{}th line processed of {}'.format(i, fname)
				tmp_block.append(line)
				if line.strip() == "M  END": 
					keyset.add(tmp_block[0].split('\t')[0])
					tmp_block = []

	data_dic = dict.fromkeys(keyset)

	for fname in list_of_parsed_files:
		with open(fname, 'rb') as r:
			tmp_block = []
			for i, line in enumerate(r):
				if i%100000==0: print '{}th line processed of {} (2nd run)'.format(i, fname)
				tmp_block.append(line)
				if line.strip() == "M  END": 
					name = tmp_block[0].strip().split('\t')
					if name[0] in zinc_dic_file:
						tmp_s = name
						tmp_s.append(zinc_dic_file[name[0]])
						data_dic[name[0]] = tmp_s[1:]
					else:
						tmp_s = name
						tmp_s.append('NO_ZINC')
						data_dic[name[0]] = tmp_s[1:]
					tmp_block = []

	with open(out_vol_cyl_id_nbs_dic, 'wb') as w:
		cPickle.dump(data_dic, w, protocol=2)

def make_zinc_dic(infile):
	keyset = set()
	prefix = '/'.join(infile.split('/')[:-1])
	with open(infile, 'rb') as r:
		for i, line in enumerate(r):
			if i%10000==0: print '{}th line processed'.format(i)
			line = line.strip().split('\t')
			keyset.add(line[0])
	with open(infile, 'rb') as r:
		zinc_dic = dict.fromkeys(keyset)
		for i, line in enumerate(r):
			if i%10000==0: print '{}th line processed (2nd run)'.format(i)
			line = line.strip().split('\t')
			zinc_dic[line[0]] = line[1]
	with open(prefix+'/ZINC_DIC_stereo.txt', 'wb') as w:
		cPickle.dump(zinc_dic, w, protocol=2)


def match_can_with_data_dic(canfile, datadic):
	keyset = set()
	keys = datadic.keys()
	with open(canfile, 'rb') as r:
		for i, line in enumerate(r):
			if i%50==0: print '{}th line processed'.format(i)
			line = line.split('\t')[0]
			keyset.add(line)
	can_dic = dict.fromkeys(keyset)
	with open(canfile, 'rb') as r:
		for i, key in enumerate(can_dic):
			if i%50==0: 
				print '{}th line processed (2nd run) {}'.format(i, key)
			try:
				can_dic[key] = datadic[key]
			except KeyError: 
				print 'KeyError {}'.format(KeyError)
	return can_dic


def select_data_from_data_dic(datadic, (min_radius, max_radius, min_height, max_height, min_volume, max_volume)):
	"""The data are stored in the dictionary as a list with the following format
	['101446140', '120.07', '3.22', '8.35', '(smi, #, A^3, A, A)', 'NO_ZINC']
	['pubchem number', volume, radius, height, units, ZINC_NUMBER]
	"""
	acc = 0
	filter_set = set()
	for keys in datadic:
		data_string = datadic[keys]
		if min_volume<=data_string[1]<=max_volume:
			if min_radius<=data_string[2]<=max_radius:
				if min_height<=data_string[3]<=max_height:
					acc+=1
					filter_set.add(keys)
	filtered_dic = dict.fromkeys(filter_set)
	for keys in filtered_dic:
		filtered_dic[keys] = datadic[keys]
		
	return acc, filtered_dic




def __main__():
	show_cylinder('c1cccc(CCCCC(C)C)c1')
	# make_zinc_dic('/home/uccahcl/Scratch/zinc15/zinc15/zinc15_merged.can')
	# zinc_dic = {}
	# with open('/home/hugues/Documents/mail/zinc/zinc/ZINC_DIC_stereo.txt', 'rb') as rb:
	# 	zinc_dic = cPickle.load(rb)


	# flist = glob.glob('/home/hugues/Documents/mail/Cylinder_processed_RKDIT/Cylinder/*part*')
	# flist = glob.glob('/home/hugues/Documents/mail/rdkit_xyz/NO_UNDESIRABLE.can_d_part_1530000.can_RD.xyz_OUT')

	# match_parsed_RDKIT_XYZ_FILE_with_ZINC(zinc_dic, flist)



	# with open('/home/hugues/Documents/mail/Cylinder_processed_RKDIT/Cylinder/DATA_DIC_STEREO.txt', 'rb') as rb:
	# 	data_dic = cPickle.load(rb)

	# can_dic = match_can_with_data_dic('/home/hugues/Documents/mail/filter/atm_nbr_split/proper_pubchem_swap/NO_UNDESIRABLE_ONLY_HETHET_ONLY_DESIRABLE.can', data_dic)
	# print len(can_dic.keys())
	# with open('/home/hugues/Documents/mail/filter/atm_nbr_split/proper_pubchem_swap/Desirable_data_dicts.can', 'rb') as r:
	# 	can_dic = cPickle.load( w)
	# select_data_from_data_dic(can_dic, (1, 2, 3, 2, 2, 3))

	# acc = 0
	# for keys in can_dic:
	# 	if acc>50:return
	# 	if can_dic[keys]==None:
	# 		acc+=1
	# 		print acc, keys, can_dic[keys]
	# 	print can_dic[keys]
__main__()
# flist = glob.glob(prefix+'*.can.xyz')
# make_smi_xyz_dictionary(flist)
# with open('/home/hugues/Documents/mail/cb-hosts/known_xyz_dict.txt', 'rb') as r:
# 	xyz_dict = pickle.load(r)

# make_cylinder_dictionary(xyz_dict, '/home/hugues/Documents/mail/cb-hosts/')

# with open('/home/hugues/Documents/mail/cb-hosts/cylinder_dict.txt', 'rb') as r:
# 	cylinder_dict = pickle.load(r)
# with open('/home/hugues/Documents/mail/cb-hosts/drawing_file.can', 'wb') as w:
# 	for k in cylinder_dict.keys():
# 		w.write('{}\t{}\t{}\t{}\n'.format(k.split('\t')[0], k.split('\t')[1],k.split('\t')[2], cylinder_dict[k][0][2]))
# np3 = merge_dicts(np1, np2)


# smi_list, tuple_list = associate_err_smi_xyz('/home/hugues/Documents/mail/cylinder_not_desirable/NO_UNDESIRABLE.can_d_part_524500.can', '/home/hugues/Documents/mail/cylinder_not_desirable/NO_UNDESIRABLE.can_d_part_524500.can.xyz')
# for items in zip(smi_list,tuple_list):
# 	print items[0], items[1].volume

# with open('/home/hugues/Documents/mail/filter/Cylinders/xyz_dict.txt', 'rb') as r:
# 	np = pickle.load(r)

# with open('/home/hugues/Documents/mail/filter/Cylinders/cylinder_dict_not_present.txt', 'rb') as r:
# 	np2 = pickle.load(r)

# keys = np2.keys()
# keys = [els.split('\t')[0] for els in keys]
# dic_clean = dict.fromkeys(keys, [])

# for key in np2:
# 	dic_clean[key.split('\t')[0]] = np2[key]

# with open('/home/hugues/Documents/mail/filter/Cylinders/cylinder_dict_not_present_smi_key.txt', 'rb') as r:
# 	np1 = pickle.load(r)

# with open('/home/hugues/Documents/mail/filter/Cylinders/cylinder_dict.txt', 'rb') as r:
# 	np2 = pickle.load(r)




# with open('/home/hugues/Documents/mail/filter/Cylinders/cylinder_file.txt', 'wb') as w:
# 	for i, key in enumerate(cyl_dic):
# 		# if i == 50: break
# 		try:
# 			w.write('{}\t{}\n'.format(key, cyl_dic[key][0][2]))
# 		except IndexError: 
# 			a+=1
# 			print '{}th, {}, bad string is {}, {}'.format(i, a, key, cyl_dic[key])

	# print cyl_dic[key][0][2]
# with open('/home/hugues/Documents/mail/filter/Cylinders/xyz_dict_tot.txt', 'wb') as w:
# 	pickle.dump(merge_dicts(np, np2), w)



# make_cylinder_dictionary(np3, prefix)



# make_smi_xyz_dictionary(flist,prefix)

# flist = [elems.split("\\")[1] for elems in flist]

