import matplotlib.pyplot as plt

def plot_amber_file(AMBERFILEOUT):
	"""
	PRE : Takes in a amber out file 
	POST : Plots it over the whole range (bar the last two entries)
	"""
	def read_amber_file(AMBERFILEOUT):
		"""
		PRE : Reads in an AMBERFILEOUT
		POST : Produces a vector to be plot 
		"""
		energies_pot = []
		with open(AMBERFILEOUT, 'rb') as r:
			for line in r:
				if 'EPtot' in line:
					energies_pot.append(line.strip().split()[-1])
		return [float(x) for x in energies_pot]

	plt.plot(read_amber_file(AMBERFILEOUT))
	plt.show()


if __name__ == "__main__":
	plot_amber_file('/media/macenrola/cb650d89-4c88-4666-a43b-08abb5756b5a/GUESTB2_prod1/outputs/GUESTB2_prod2/CB_prod2.out')