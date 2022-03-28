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
		x += coord[0]*mass
		y += coord[1]*mass
		z += coord[2]*mass
		totmass += mass
	global cM
	cM = numpy.array([x/totmass, y/totmass, z/totmass])


	I = []
	for index in range(9):
		I.append(0)

	for (mass, coord) in zip(masses, coords):
		temp_x, temp_y, temp_z = coord[0], coord[1], coord[2]
		temp_x -= x
		temp_y -= y
		temp_z -= z

		I[0] += mass * (temp_y**2 + temp_z**2)
		I[1] -= mass * temp_x * temp_y
		I[2] -= mass * temp_x * temp_z
		I[3] -= mass * temp_x * temp_y
		I[4] += mass * (temp_x**2 + temp_z**2)
		I[5] -= mass * temp_y * temp_z
		I[6] -= mass * temp_x * temp_z
		I[7] -= mass * temp_y * temp_z
		I[8] += mass * (temp_x**2 + temp_y**2)

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

def main():
	align_protein_to_inertial_axes("prot.pdb", "prot_al.pdb")
	([minX, minY, minZ],[maxX, maxY, maxZ]) = get_protein_dimensions("prot_al.pdb")
	print("The Inertia Axis Aligned Bounding Box (IABB) dimensions are (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))
	print("The Inertia Axis Aligned Bounding Box (IABB) volume is %.2f A3" % ((maxX-minX)*(maxY-minY)*(maxZ-minZ)))

if __name__ == "__main__":
    main()



