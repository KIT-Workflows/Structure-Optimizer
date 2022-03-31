import sys, yaml, os, tarfile

from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints, rdmolfiles
from rdkit.ML.Cluster import Butina


import openbabel
from openbabel import pybel

def gen_conformers(mol, numConfs=100, maxAttempts=1000, pruneRmsThresh=0.1, useExpTorsionAnglePrefs=True, useBasicKnowledge=True, enforceChirality=True):
	ids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, maxAttempts=maxAttempts, pruneRmsThresh=pruneRmsThresh, useExpTorsionAnglePrefs=useExpTorsionAnglePrefs, useBasicKnowledge=useBasicKnowledge, enforceChirality=enforceChirality, numThreads=8)
	return list(ids)
	
def write_conformers_to_sdf(mol, filename, rmsClusters, conformerPropsDict, minEnergy):
	w = Chem.SDWriter(filename)
	for cluster in rmsClusters:
		for confId in cluster:
			for name in mol.GetPropNames():
				mol.ClearProp(name)
			conformerProps = conformerPropsDict[confId]
			mol.SetIntProp("conformer_id", confId + 1)
			for key in conformerProps.keys():
				mol.SetProp(key, str(conformerProps[key]))
			e = conformerProps["energy_abs"]
			if e:
				mol.SetDoubleProp("energy_delta", e - minEnergy)
			w.write(mol, confId=confId)
	w.flush()
	w.close()
	
def calc_energy(mol, conformerId, minimizeIts):
	ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conformerId)
	ff.Initialize()
	ff.CalcEnergy()
	results = {}
	if minimizeIts > 0:
		results["converged"] = ff.Minimize(maxIts=minimizeIts)
	results["energy_abs"] = ff.CalcEnergy()
	return results
	
def cluster_conformers(mol, mode="RMSD", threshold=2.0):
	if mode == "TFD":
		dmat = TorsionFingerprints.GetTFDMatrix(mol)
	else:
		dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
	rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
	return rms_clusters
	
def align_conformers(mol, clust_ids):
	rmslist = []
	AllChem.AlignMolConformers(mol, confIds=clust_ids, RMSlist=rmslist)
	return rmslist
		

if __name__ == '__main__':

	with open('rendered_wano.yml') as file:
		wano_file = yaml.full_load(file)

	input_file = wano_file["Input-File"]
	# if len(sys.argv) < 4:
	# 	print ("Usage: conf_gen.py <sdf input> <num conformers> <max attempts> <prune threshold> <cluster method: (RMSD|TFD) = RMSD> <cluster threshold = 0.2> <minimize iterations: = 0>")
	# 	exit()
	# #python test1.py structure.sdf 40 10 0.5
	# input_file = sys.argv[1]

	mol = pybel.readfile("xyz", input_file)
	os.system("obabel " + input_file + " -O structure.sdf")	
	input_file = "structure.sdf"
	
	numConfs = int(wano_file["Parameters"]["nconf"]) #int(sys.argv[2])
	maxAttempts = int(wano_file["Parameters"]["max-attempts"])#int(sys.argv[3])
	pruneRmsThresh = wano_file["Parameters"]["cluster-thr"] #float(sys.argv[4])
	clusterMethod = wano_file["Parameters"]["cluster-algo"]
	clusterThreshold = wano_file["RMSD-thr"]
	minimizeIterations = wano_file["Parameters"]["max-opt"] 

	suppl = Chem.ForwardSDMolSupplier(input_file)
	i=0
	for mol in suppl:
		i = i+1
		if mol is None: continue
		m = Chem.AddHs(mol)
		# generate the confomers
		conformerIds = gen_conformers(m, numConfs, maxAttempts, pruneRmsThresh, True, True, True)
		conformerPropsDict = {}
		for conformerId in conformerIds:
			# energy minimise (optional) and energy calculation
			props = calc_energy(m, conformerId, minimizeIterations)
			conformerPropsDict[conformerId] = props
		# cluster the conformers
		rmsClusters = cluster_conformers(m, clusterMethod, clusterThreshold)

		print ("Molecule", i, ": generated", len(conformerIds), "conformers and", len(rmsClusters), "clusters")
		rmsClustersPerCluster = []
		clusterNumber = 0
		minEnergy = 9999999999999
		for cluster in rmsClusters:
			clusterNumber = clusterNumber+1
			rmsWithinCluster = align_conformers(m, cluster)
			for conformerId in cluster:
				e = props["energy_abs"]
				if e < minEnergy:
					minEnergy = e
				props = conformerPropsDict[conformerId]
				props["cluster_no"] = clusterNumber
				props["cluster_centroid"] = cluster[0] + 1
				idx = cluster.index(conformerId)
				if idx > 0:
					props["rms_to_centroid"] = rmsWithinCluster[idx-1]
				else:
					props["rms_to_centroid"] = 0.0

		write_conformers_to_sdf(m, "my_conf_" + str(i) + ".sdf", rmsClusters, conformerPropsDict, minEnergy)

	i=0
	for mol in pybel.readfile("sdf", "my_conf_1.sdf"):
		mol.write("xyz", "mol_"+str(i)+".xyz")
		i+=1

	archive = tarfile.open("mol.tar.gz", "w|gz")

	for fname in os.listdir():
		if fname.startswith("mol_"):
			archive.add(fname, arcname=fname)
			os.remove(os.path.join(fname)) 

	archive.close()

	#os.system("obabel *.sdf -oxyz -m")
	#os.system('obabel my_conf_1.sdf -O output.xyz --split')
	# struct = cmd.load('my_conf_1.sdf')
	# cmd.save('my_file.xyz')
	# new_struct = rdmolfiles.MolToXYZFile(struct)


	# from pymatgen.core import Molecule
	# #from pymatgen.io.babel import confab_conformers
	# #from openbabel import openbabel
	# from openbabel import openbabel as ob
	# from openbabel import pybel as pb

	# mol_struct = Molecule.from_file("Henrik6.xyz")
	# pb.Molecule.conformers(forcefield='mmff94', freeze_atoms=None, rmsd_cutoff=0.5, energy_cutoff=50.0, conf_cutoff=100000, verbose=False)


	# #print(mol_struct)

