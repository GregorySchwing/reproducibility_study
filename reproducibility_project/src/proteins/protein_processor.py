from __future__ import print_function

__author__  = 'Pablo Guardado Calvo'
__version__ = '0.1'
__email__   = 'pablo.guardado (at) gmail.com'
__date__    = '13/08/2015'


"""Methods used to create process pdb file job statepoint."""


"""
Calculate and display the dimensions of a protein.


This is a first version, please use at your own risk!

REQUIREMENTS

numpy (http://numpy.scipy.org) that should be built into the newers versions of Pymol

(c) Pablo Guardado Calvo

Based on "inertia_tensor.py" (c) 2010 by Mateusz Maciejewski

License: MIT

"""

# For a dimer, Myc-Max, aligning decreases the box volume by 85721.536846008 A^3

from prody.proteins.pdbfile import parsePDB
from prody.proteins.pdbfile import writePDB
from prody.measure.transform import *
from prody.measure.measure import *

import numpy

def matriz_inercia(selection):
	'''
	DESCRIPTION

	The method calculates the mass center, the inertia tensor and the eigenvalues and eigenvectors
	for a given selection. Mostly taken from inertia_tensor.py
	'''

	masses = selection.getMasses()
	coords = selection.getCoords()
	totmass = 0.0
	x,y,z = 0,0,0
	for (mass, coord) in zip(masses, coords):
		m = mass
		m = 1
		x += coord[0]*m
		y += coord[1]*m
		z += coord[2]*m
		totmass += mass
	global cM
	cM = numpy.array([x/totmass, y/totmass, z/totmass])


	I = []
	for index in range(9):
		I.append(0)

	for (mass, coord) in zip(masses, coords):
		m = mass
		m = 1
		temp_x, temp_y, temp_z = coord[0], coord[1], coord[2]
		temp_x -= x
		temp_y -= y
		temp_z -= z

		I[0] += m * (temp_y**2 + temp_z**2)
		I[1] -= m * temp_x * temp_y
		I[2] -= m * temp_x * temp_z
		I[3] -= m * temp_x * temp_y
		I[4] += m * (temp_x**2 + temp_z**2)
		I[5] -= m * temp_y * temp_z
		I[6] -= m * temp_x * temp_z
		I[7] -= m * temp_y * temp_z
		I[8] += m * (temp_x**2 + temp_y**2)

	global tensor
	tensor = numpy.array([(I[0:3]), (I[3:6]), (I[6:9])])

	global autoval, autovect, ord_autoval, ord_autovect
	autoval, autovect = numpy.linalg.eig(tensor)
	auto_ord = numpy.argsort(autoval)
	ord_autoval = autoval[auto_ord]
	ord_autovect_complete = autovect[:, auto_ord].T
	ord_autovect = numpy.around(ord_autovect_complete, 3)

	return ord_autoval



def translacion_cM(selection):
	'''
	DESCRIPTION

	Translate the center of mass of the molecule to the origin.
	'''
	masses = selection.getMasses()
	coords = selection.getCoords()
	totmass = 0.0
	x,y,z = 0,0,0
	for (mass, coord) in zip(masses, coords):
		m = mass
		m = 1
		x += coord[0]*m
		y += coord[1]*m
		z += coord[2]*m
		totmass += m
	cM = numpy.array([x/totmass, y/totmass, z/totmass])
	trans_array = ([1, 0, 0, -cM[0], 0, 1, 0, -cM[1], 0, 0, 1, -cM[2], 0, 0, 0, 1])
	trans_array_2D = numpy.reshape(trans_array, (-1, 4))
	model_trans = Transformation(trans_array_2D)
	newMol = applyTransformation(model_trans, selection)

def align_protein_to_inertial_axes(
	infilename,
	outfilename,
):
	molecule = parsePDB(infilename)
	selection = molecule.protein
	translacion_cM(selection)
	matriz_inercia(selection)
	global transf, transf_array, ord_autovect_array, transf_array_print
	ord_autovect_array = numpy.array([[ord_autovect[0][0], ord_autovect[0][1], ord_autovect[0][2]],
	                                  [ord_autovect[1][0], ord_autovect[1][1], ord_autovect[1][2]],
					  [ord_autovect[2][0], ord_autovect[2][1], ord_autovect[2][2]]])
	if numpy.linalg.det(ord_autovect_array) == -1:
		ord_autovect_array = numpy.array([[ord_autovect[2][0], ord_autovect[2][1], ord_autovect[2][2]],
	                                  	  [ord_autovect[1][0], ord_autovect[1][1], ord_autovect[1][2]],
	 				  	  [ord_autovect[0][0], ord_autovect[0][1], ord_autovect[0][2]]])
	transf = numpy.transpose(ord_autovect_array)
	transf_array = numpy.array([transf[0][0], transf[0][1], transf[0][2], 0,
	                            transf[1][0], transf[1][1], transf[1][2], 0,
				    transf[2][0], transf[2][1], transf[2][2], 0,
				    0, 0, 0, 1])
	transf_array_2D = numpy.reshape(transf_array, (-1, 4))
	model_transf = Transformation(transf_array_2D)
	newMol = applyTransformation(model_transf, selection)
	molecule = writePDB(outfilename,selection)

def get_protein_dimensions(
    infilename,
):	
	molecule = parsePDB(infilename)
	selection = molecule.protein
	coords = selection.getCoords()
	x = coords[:,0]
	y = coords[:,1]
	z = coords[:,2]
	minX = numpy.min(x)
	minY = numpy.min(y)
	minZ = numpy.min(z)
	maxX = numpy.max(x)
	maxY = numpy.max(y)
	maxZ = numpy.max(z)
	return ([minX, minY, minZ],[maxX, maxY, maxZ])

def get_protein_path(
    infilename,
):	
	from src import proteins
	from os.path import join
	from os.path import abspath
	from os.path import dirname
	prot_path = (
		join(str(dirname(abspath(proteins.__file__))), infilename)
	)
	return prot_path

def get_total_charge(
    job,
):
	from prody.trajectory.psffile import parsePSF, writePSF
	from prody.proteins.pdbfile import parsePDB, writePDB
	getCharges

def merge_solv_and_solute(
    job,
):	

	from numpy import zeros
	from prody.trajectory.psffile import parsePSF, writePSF
	from prody.proteins.pdbfile import parsePDB, writePDB
	from prody.measure.transform import moveAtoms
	from parmed.charmm import CharmmPsfFile
	from parmed.formats.registry import load_file
	import numpy as np
	import pandas as pd
	#Shuift geometric center of solvent box to origin
	print("Centering solvent on origin")
	uncentered_positions = parsePDB(job.fn("mosdef_box_0.pdb"))
	moveAtoms(uncentered_positions, to=zeros(3), ag=True)
	writePDB(job.fn("mosdef_box_0_centered.pdb"), uncentered_positions)

	#Load centered sovent box into Partmed
	solvbox2 = CharmmPsfFile(job.fn("mosdef_box_0.psf"))
	solvbox2_pos = load_file(job.fn("mosdef_box_0_centered.pdb"))


	#Load protein structure into Partmed
	pro2_parm = CharmmPsfFile(get_protein_path(job.sp.pdbid+".psf"))

	pro2_parm_pos = load_file(get_protein_path(job.sp.pdbid+"_aligned.pdb"))
	parmPosComb = pro2_parm_pos + solvbox2_pos

	print("Superimposing protein and solvent")
	#Combine topolgies in parmed and write intermediate files
	# Also note that the protein has to be on the left when merging (prot+solv)
	# Very important not ro renumber since we need to cross reference with the psf.
	systemPSF = CharmmPsfFile.from_structure(pro2_parm + solvbox2)
	systemPSF.write_psf(job.fn("intermediate.psf"))
	parmPosComb.write_pdb(job.fn("intermediate.pdb"), renumber=False, use_hetatoms=False)

	#Load intermediate structure into Partmed
	systemPSFComb = CharmmPsfFile(job.fn("intermediate.psf"))
	systemPSFComb_pos = load_file(job.fn("intermediate.pdb"))

	df = systemPSFComb_pos.to_dataframe()
	df2 = systemPSFComb.to_dataframe()

	# Separate protein from solvent using Mosdef's SYS segment ID
	# Very important not to change this before merging topologies
	endOfProtein = df2.segid.eq('SYS').idxmax()
	startOfSolvent = endOfProtein+1
	endOfSolvent = len(df2.index)
	prot = "@0-{end}".format(end=endOfProtein)
	solv = "@{start}-{end}".format(start=startOfSolvent, end=endOfSolvent)

	print("Masking solvent residues possessing an atom within 2.4 A of protein")
	amberMaskRes = "{protS}<:2.4&{solvS}".format(protS=prot, solvS=solv)
	badWatersRes = systemPSFComb_pos[amberMaskRes]

	# Create dataframes for checking existence unique combinations in bad water subset. 
	bwDF = badWatersRes.to_dataframe()
	allAtoms = df2[["name","resname","resnum"]]
	badWaters = bwDF[["name","resname","resnum"]]

	# Create a list of length atoms where values are true and false, indicating whether this atom should be removed.
	df = pd.merge(allAtoms, badWaters, on=["name","resname","resnum"], how='left', indicator='Exist')
	df['Exist'] = np.where(df.Exist == 'both', True, False)
	stripInput = df['Exist'].to_numpy()

	print("Stripping masked solvent")
	# Strip bad waters from coordinates object
	splitRes = systemPSFComb.strip(stripInput)
	systemPSFComb.write_psf(job.fn("solvated.psf"))

	print("Writing combined topology/coordinate files")
	# Strip bad waters from topology object
	splitRes = systemPSFComb_pos.strip(stripInput)
	systemPSFComb_pos.write_pdb(job.fn("solvated.pdb"), use_hetatoms=False)

# Convert segment SYS to W00 to WNN, so prevent repeated (resid, segment) pairs
def fix_segment(
    job,
):	
	Text_File_Import = open(job.fn("solvated.psf"), 'r')

	Text_lines = Text_File_Import.readlines()
	counter = 00
	segID = "W00"
	lastResid = 0
	with open(job.fn("solvated.psf"), 'w') as the_file:	
		for line in Text_lines:
			User_Inputs = line.split()
			resName = 0
			try:
				resName = User_Inputs.index("WAT")
			except ValueError:
				"print no water"		

			if (resName > 0):
				if (User_Inputs[resName-1] == str(1) and lastResid == str(9999)):
					counter = counter + 1
					segID = "W{0:0=2d}".format(counter)
				newLine = line.replace('SYS', segID)
				the_file.write(newLine)
				lastResid = User_Inputs[resName-1]
			else:
				the_file.write(line)

def compute_ion_numbers_and_positions(job):

	from vmd import evaltcl
	template = get_protein_path("ionizeTemplate.tcl")
	helperMethods = get_protein_path("concentrationHelpers.tcl")
	# Read in the file
	with open(template, 'r') as file :
		filedata = file.read()

	# Replace the target string
	filedata = filedata.replace("HELPER_METHODS", helperMethods)
	filedata = filedata.replace("PATH_2_PROT_PSF", get_protein_path(job.sp.pdbid+".psf"))
	filedata = filedata.replace("PATH_2_PROT_PDB", get_protein_path(job.sp.pdbid+"_aligned.pdb"))
	filedata = filedata.replace("PATH_2_SOLVATED_PSF", job.doc.path_to_solvated_psf)
	filedata = filedata.replace("PATH_2_SOLVATED_PDB", job.doc.path_to_solvated_pdb)
	filedata = filedata.replace("CATION_NAME", job.sp.cat_name)
	filedata = filedata.replace("CATION_VAL", str(job.sp.cat_val))
	filedata = filedata.replace("ANION_NAME", job.sp.an_name)
	filedata = filedata.replace("ANION_VAL", str(job.sp.an_val))
	filedata = filedata.replace("SALC_CONC_VALUE", str(job.sp.salt_conc))
	# Write the file out again
	with open(job.fn("filled_template.tcl"), 'w') as file:
		file.write(filedata)


	print("Computing number and positions of ions in", job)
	ions = evaltcl("source " + job.fn("filled_template.tcl"))
	ionsList = ions.split()

	job.document["NCATION"] = ionsList[0]
	job.document["NANION"] = ionsList[1]
	with open(job.fn("NUMIONS.txt"), "w") as file:
		file.write(ions + "\n")

def build_ions_psf(job):
	import mbuild as mb
	from foyer import Forcefield
	import mbuild.formats.charmm_writer as mf_charmm
	import mbuild.formats.gomc_conf_writer as gomc_control
	from reproducibility_project.src.utils.forcefields import get_ff_path_ion

	liq_box = mb.Box([1.0, 1.0, 1.0])

	ions_2_smiles = {
		"SOD": '[Na]',
		"CLA": '[Cl-]'
	}

	cation = mb.load(ions_2_smiles[job.sp.cat_name], smiles=True)
	cation.name = job.sp.cat_name
	num_cations = int(job.document["NCATION"])
	FF_file_cation = get_ff_path_ion("custom", job.sp.cat_name)

	filled_cation_box = mb.fill_box(
		compound=[cation],
		n_compounds=[num_cations],
		box=liq_box
	)

	cations_list = [cation.name]

	# Build charmm object for cations
	charmm_cations = mf_charmm.Charmm(filled_cation_box,
							job.sp.cat_name,
							structure_box_1=None,
							filename_box_1=None,
							ff_filename=job.sp.cat_name,
							forcefield_selection=FF_file_cation,
							residues=cations_list,
							bead_to_atom_name_dict=None,
							fix_residue=None,
							gomc_fix_bonds_angles=None,
							reorder_res_in_pdb_psf=True
							)

	charmm_cations.write_inp()
	charmm_cations.write_psf()


	anion = mb.load(ions_2_smiles[job.sp.an_name], smiles=True)
	anion.name = job.sp.an_name
	num_anions = int(job.document["NANION"])
	FF_file_anion = get_ff_path_ion("custom", job.sp.an_name)

	filled_anion_box = mb.fill_box(
		compound=[anion],
		n_compounds=[num_anions],
		box=liq_box
	)

	anions_list = [anion.name]


	# Build charmm object for anions
	charmm_anions = mf_charmm.Charmm(filled_anion_box,
							job.sp.an_name,
							structure_box_1=None,
							filename_box_1=None,
							ff_filename=job.sp.an_name,
							forcefield_selection=FF_file_anion,
							residues=anions_list,
							bead_to_atom_name_dict=None,
							fix_residue=None,
							gomc_fix_bonds_angles=None,
							reorder_res_in_pdb_psf=True
							)

	charmm_anions.write_inp()
	charmm_anions.write_psf()


def merge_ions_and_system(
    job,
):	

	from numpy import zeros
	from prody.trajectory.psffile import parsePSF, writePSF
	from prody.proteins.pdbfile import parsePDB, writePDB
	from prody.measure.transform import moveAtoms
	from parmed.charmm import CharmmPsfFile
	from parmed.formats.registry import load_file
	import numpy as np
	import pandas as pd
	import re
	import os  

	p = re.compile("(\w+?)-ions_(\d+)-(\w+?).pdb")

	from os import listdir
	ionFiles = [f for f in listdir(job.ws) if ((job.isfile(f) and p.match(f)))]

	print(ionFiles)
	for ionFile in ionFiles:
		if (job.sp.an_name in ionFile):
			os.rename(ionFile,job.sp.an_name+'.pdb')
		elif (job.sp.cat_name in ionFile):
			os.rename(ionFile,job.sp.cat_name+'.pdb')

	#Load centered sovent box into Partmed
	print("Loading system")
	system = CharmmPsfFile(job.doc.path_to_solvated_psf)
	system_pos = load_file(job.doc.path_to_solvated_pdb)
	parmPosComb = system_pos

	#Shuift geometric center of solvent box to origin
	print("Merging ions with system")
	for f in ionFiles:
		ions = load_file(job.fn(f))
		parmPosComb = parmPosComb + ions
		system = system + ions

	parmPosComb.write_pdb(job.fn("combined.pdb"), renumber=False, use_hetatoms=False)
	system.write_psf(job.fn("combined.psf"))
	quit()
	#Load centered sovent box into Partmed
	solvbox2 = CharmmPsfFile(job.fn("mosdef_box_0.psf"))
	solvbox2_pos = load_file(job.fn("mosdef_box_0_centered.pdb"))


	#Load protein structure into Partmed
	pro2_parm = CharmmPsfFile(get_protein_path(job.sp.pdbid+".psf"))

	pro2_parm_pos = load_file(get_protein_path(job.sp.pdbid+"_aligned.pdb"))
	parmPosComb = pro2_parm_pos + solvbox2_pos

	print("Superimposing protein and solvent")
	#Combine topolgies in parmed and write intermediate files
	# Also note that the protein has to be on the left when merging (prot+solv)
	# Very important not ro renumber since we need to cross reference with the psf.
	systemPSF = CharmmPsfFile.from_structure(pro2_parm + solvbox2)
	systemPSF.write_psf(job.fn("intermediate.psf"))
	parmPosComb.write_pdb(job.fn("intermediate.pdb"), renumber=False, use_hetatoms=False)

	#Load intermediate structure into Partmed
	systemPSFComb = CharmmPsfFile(job.fn("intermediate.psf"))
	systemPSFComb_pos = load_file(job.fn("intermediate.pdb"))

	df = systemPSFComb_pos.to_dataframe()
	df2 = systemPSFComb.to_dataframe()

	# Separate protein from solvent using Mosdef's SYS segment ID
	# Very important not to change this before merging topologies
	endOfProtein = df2.segid.eq('SYS').idxmax()
	startOfSolvent = endOfProtein+1
	endOfSolvent = len(df2.index)
	prot = "@0-{end}".format(end=endOfProtein)
	solv = "@{start}-{end}".format(start=startOfSolvent, end=endOfSolvent)

	print("Masking solvent residues possessing an atom within 2.4 A of protein")
	amberMaskRes = "{protS}<:2.4&{solvS}".format(protS=prot, solvS=solv)
	badWatersRes = systemPSFComb_pos[amberMaskRes]

	# Create dataframes for checking existence unique combinations in bad water subset. 
	bwDF = badWatersRes.to_dataframe()
	allAtoms = df2[["name","resname","resnum"]]
	badWaters = bwDF[["name","resname","resnum"]]

	# Create a list of length atoms where values are true and false, indicating whether this atom should be removed.
	df = pd.merge(allAtoms, badWaters, on=["name","resname","resnum"], how='left', indicator='Exist')
	df['Exist'] = np.where(df.Exist == 'both', True, False)
	stripInput = df['Exist'].to_numpy()

	print("Stripping masked solvent")
	# Strip bad waters from coordinates object
	splitRes = systemPSFComb.strip(stripInput)
	systemPSFComb.write_psf(job.fn("solvated.psf"))

	print("Writing combined topology/coordinate files")
	# Strip bad waters from topology object
	splitRes = systemPSFComb_pos.strip(stripInput)
	systemPSFComb_pos.write_pdb(job.fn("solvated.pdb"), use_hetatoms=False)

def ionize(
    job,
):	
	compute_ion_numbers_and_positions(job)
	build_ions_psf(job)
	merge_ions_and_system(job)


def main():
	align_protein_to_inertial_axes("prot.pdb", "prot_al.pdb")
	([minX, minY, minZ],[maxX, maxY, maxZ]) = get_protein_dimensions("prot_al.pdb")
	print("The Inertia Axis Aligned Bounding Box (IABB) dimensions are (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))
	print("The Inertia Axis Aligned Bounding Box (IABB) volume is %.2f A3" % ((maxX-minX)*(maxY-minY)*(maxZ-minZ)))


if __name__ == "__main__":
    main()



