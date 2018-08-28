import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
from subprocess import call
sys.path.append('/home/macenrola/Thesis/Molecule_operations/')
from mol_ops import make_pdb_complex_with_named_residues, align_mol

def one_generate_pdb_script(name, fout=''):
	"""
	PRE:  Takes in a pqr mol
	POST: Generates an apbs script that is meant to obtain both its polar and no polar contribution to binding energy
	"""
	name = name.split('/')[-1]
	if fout == '':
		fout = name+'_solv_script'
	head_script = [
	'read',
	'  mol pqr {}'.format(name),
	'end'
	]
	bodyscript=[
	'elec name {}-solv # Electrostatics calculation on the solvated state'.format(name),
	'  mg-manual # Specify the mode for APBS to run',
	'  dime 97 97 97 # The grid dimensions',
	'  nlev 4 # Multigrid level parameter',
	'  grid 0.33 0.33 0.33 # Grid spacing',
	'  gcent mol {} # Center the grid on molecule 1'.format(1),
	'  mol {} # Perform the calculation on molecule 1'.format(1),
	'  lpbe # Solve the linearized Poisson-Boltzmann equation',
	'  bcfl mdh # Use all multipole moments when calculating the potential',
	'  pdie 2.25 # Solute dielectric',
	'  sdie 78.54 # Solvent dielectric',
	'  chgm spl2 # Spline-based discretization of the delta functions',
	'  srfm mol # Molecular surface definition',
	'  srad 1.4 # Solvent probe radius (for molecular surface)',
	'  swin 0.3 # Solvent surface spline window (not used here)',
	'  sdens 10.0 # Sphere density for accessibility object',
	'  temp 298.15 # Temperature',
	'  calcenergy total # Calculate energies',
	'  calcforce no # Do not calculate forces',
	'end',
	'elec name {}-ref # Calculate potential for reference (vacuum) state'.format(name),
	'  mg-manual',
	'  dime 97 97 97',
	'  nlev 4',
	'  grid 0.33 0.33 0.33',
	'  gcent mol {}'.format(1),
	'  mol {}'.format(1),
	'  lpbe',
	'  bcfl mdh',
	'  pdie 2.25',
	'  sdie 1.0',
	'  chgm spl2',
	'  srfm mol',
	'  srad 1.4',
	'  swin 0.3',
	'  sdens 10.0',
	'  temp 298.15',
	'  calcenergy total',
	'  calcforce no',
	'end',

	'apolar name {}-apol'.format(name),
	'  bconc 0',
	'  gamma 0.0418',
	'  press 0',
	'  grid 0.33 0.33 0.33 # Grid spacing',
	'  mol {}'.format(1),
	'  sdens 10.0',
	'  srad 1.4',
	'  srfm sacc',
	'  swin 0.3',
	'  calcenergy total',
	'  calcforce no',
	'  dpos 0.2',
	'  temp 300.0',
	'end'
	]
	head_script.extend(bodyscript)

	bottom_script = [
	'print ElecEnergy {0}-solv - {0}-ref end'.format(name),
	'print ApolEnergy {}-apol end'.format(name),
	]
	head_script.extend(bottom_script)
	head_script.append('print ElecEnergy {0}-solv - {0}-ref end'.format(name))
	head_script.append('print ApolEnergy {0}-apol end'.format(name))
	with open(fout, 'wb') as w:
		w.write('\n'.join(head_script))
	return 

def generate_pdb_script(guestmol, name='GUEST.pqr', fout='solvation_apbs'):
	"""
	PRE:  Takes in a rdkit mol
	POST: Generates an apbs script that is meant to obtain both its polar and no polar contribution to binding energy
	"""

	mols = {1:name, 2:'CB.pqr', 3:'CB-{}'.format(name)}
	head_script = [
	'read',
	'  mol pqr {}'.format(name),
	'  mol pqr {}'.format('CB.pqr'),
	'  mol pqr {}'.format('CB-{}'.format(name)),
	'end'
	]
	for i in range(3):
		bodyscript=[
		'elec name {}-solv # Electrostatics calculation on the solvated state'.format(mols[i+1]),
		'  mg-manual # Specify the mode for APBS to run',
		'  dime 97 97 97 # The grid dimensions',
		'  nlev 4 # Multigrid level parameter',
		'  grid 0.33 0.33 0.33 # Grid spacing',
		'  gcent mol {} # Center the grid on molecule 1'.format(i+1),
		'  mol {} # Perform the calculation on molecule 1'.format(i+1),
		'  lpbe # Solve the linearized Poisson-Boltzmann equation',
		'  bcfl mdh # Use all multipole moments when calculating the potential',
		'  pdie 2.25 # Solute dielectric',
		'  sdie 78.54 # Solvent dielectric',
		'  chgm spl2 # Spline-based discretization of the delta functions',
		'  srfm mol # Molecular surface definition',
		'  srad 1.4 # Solvent probe radius (for molecular surface)',
		'  swin 0.3 # Solvent surface spline window (not used here)',
		'  sdens 10.0 # Sphere density for accessibility object',
		'  temp 298.15 # Temperature',
		'  calcenergy total # Calculate energies',
		'  calcforce no # Do not calculate forces',
		'end',
		'elec name {}-ref # Calculate potential for reference (vacuum) state'.format(mols[i+1]),
		'  mg-manual',
		'  dime 97 97 97',
		'  nlev 4',
		'  grid 0.33 0.33 0.33',
		'  gcent mol {}'.format(i+1),
		'  mol {}'.format(i+1),
		'  lpbe',
		'  bcfl mdh',
		'  pdie 2.25',
		'  sdie 1.0',
		'  chgm spl2',
		'  srfm mol',
		'  srad 1.4',
		'  swin 0.3',
		'  sdens 10.0',
		'  temp 298.15',
		'  calcenergy total',
		'  calcforce no',
		'end',

		'apolar name {}-apol'.format(mols[i+1]),
		'  bconc 0',
		'  gamma 0.0418',
		'  press 0',
		'  grid 0.33 0.33 0.33 # Grid spacing',
		'  mol {}'.format(i+1),
		'  sdens 10.0',
		'  srad 1.4',
		'  srfm sacc',
		'  swin 0.3',
		'  calcenergy total',
		'  calcforce no',
		'  dpos 0.2',
		'  temp 300.0',
		'end'
		]
		head_script.extend(bodyscript)
	for i in range(3):
		bottom_script = [
		'print ElecEnergy {0}-solv - {0}-ref end'.format(mols[i+1]),
		'print ApolEnergy {}-apol end'.format(mols[i+1]),
		]
		head_script.extend(bottom_script)
	head_script.append('print ElecEnergy {0}-solv - {0}-ref - {1}-solv + {1}-ref - {2}-solv + {2}-ref end'.format(mols[3], mols[2], mols[1]))
	head_script.append('print ApolEnergy {0}-apol - {1}-apol - {2}-apol end'.format(mols[3], mols[2], mols[1]))
	with open(fout, 'wb') as w:
		w.write('\n'.join(head_script))

	# print '\n'.join(head_script)
	return 

def convert_pdb_to_pqr(guestpathpdb, cbpathpdb, complexpathpdb):
	"""
	PRE  : Takes the path to the guest.PDB, cb.PDB and complex.PDB
		   Requires Openbabel 2.3.2 or later 
	POST : Generates the PQR files along them by replacing the ".PDB" by ".PQR"
	"""
	for names in [guestpathpdb, cbpathpdb, complexpathpdb]:
		cmd =  '{} {} -O {}'.format(pdb2pqr_cmd, names, names.replace('.pdb', '.pqr'))
		print cmd
		call(cmd, shell=True)

def one_convert_pdb_to_pqr(molpdb, 	pdb2pqr_cmd = '/usr/bin/obabel'):
	"""
	PRE  : Takes the path to the mol.PDB
		   Requires Openbabel 2.3.2 or later 
	POST : Generates the PQR files along them by replacing the ".PDB" by ".PQR"
	"""
	molpqr = molpdb[:-4]+'.pqr'
	cmd =  '{} {} -O {}'.format(pdb2pqr_cmd, molpdb, molpqr)
	print cmd
	call(cmd, shell=True)
	return

def get_solvation_for_pdb(molpdbfile,apbs_cmd = '/usr/bin/apbs'):
	"""
	PRE  : Takes in a file in pdb format containing a molecule
		   Requires APBS
	POST : Computes the solvation polar and apolar of the molecule and returns the numbers
		   Generates a proper APBS script and calls it, after having converted the molecule to a proper pqr file
		   Returns the free energies in kcal/mol as a tuple, (APOL, ELEC)
	"""
	molpqr = molpdbfile[:-4]+'.pqr'
	apbs_script = molpdbfile[:-4] + '_solv_script'
	apbs_output = apbs_script + '_output'
	one_convert_pdb_to_pqr(molpdbfile)
	one_generate_pdb_script(molpqr,fout=apbs_script)
	call('cd {3} && {0} {1} > {2}'.format(apbs_cmd, apbs_script.split('/')[-1], apbs_output.split('/')[-1], '/'.join(apbs_output.split('/')[:-1])), shell=True)
	with open(apbs_output) as r:
		lines = r.readlines()
	for s in lines[::-1]:
		if '  Global net APOL energy' in s:
			apol_s = s
		elif '  Global net ELEC energy' in s:
			elec_s = s
			break
	apol_energy = apol_s.split()[-2]
	elec_energy = elec_s.split()[-2]

	print 'solvation energy (APOL+ELEC) is {}+{}={} kJ/mol (or {}+{}={} kcal/mol)'.format(float(apol_energy), float(elec_energy), float(apol_energy)+float(elec_energy),float(apol_energy)/4.184, float(elec_energy)/4.184, (float(apol_energy)+float(elec_energy))/4.184)
	return float(apol_energy)/4.184, float(elec_energy)/4.184

if __name__=="__main__":
	pass
	get_solvation_for_pdb('/home/macenrola/Thesis/hydrocarbons/CB_REFERENCE_VALUES/CB_candidate.pdb')
	# pdb2pqr_cmd = '/usr/bin/pdb2pqr'
	# pdb2pqr_cmd = '/usr/bin/obabel'
	# apbs_cmd = '/usr/bin/apbs'

	# #############################
	# # Solve the charge issue
	# #############################

	# print get_solvation_for_pdb('/home/macenrola/Thesis/hydrocarbons/splitby1hydrocarbon-19/CB_candidate-SDF2PDB.pdb')
	# smi = 'C[N+](C12CC3C4C(C2)C2C(C1)C3CC(C4)(C2)[N+](C)(C)C)(C)C'
	# smi = '[NH4+]'
	# smi = 'C1=CC=NC=C1'
	# smi = 'c1ccccc1'
	# smi = 'C1=CC(=CC(=C1)O)O'
	# # smi = 'CC(=O)C'
	# mol = Chem.MolFromSmiles(smi)
	# mol = Chem.AddHs(mol)
	# AllChem.EmbedMolecule(mol, AllChem.ETKDG())
	# mol = align_mol(Chem.MolToMolBlock(mol))
	# make_pdb_complex_with_named_residues(mol, './GUEST.pdb', './CB.pdb', './CB-GUEST.pdb', isCB=False)
	# convert_pdb_to_pqr('./GUEST.pdb', './CB.pdb', './CB-GUEST.pdb')
	# generate_pdb_script(mol, 'GUEST.pqr')
	# call_apbs('solvation_apbs')
	# get_solvation_for_pdb('')