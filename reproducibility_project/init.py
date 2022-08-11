"""Initialize signac statepoints."""
import itertools
import os

import numpy as np
import signac
import unyt as u
from numpy import ModuleDeprecationWarning
from numpy import ModuleDeprecationWarning

from src.proteins.protein_processor import (
    align_protein_to_inertial_axes,
    get_protein_dimensions,
    get_protein_path,
)

def dict_product(dd):
    """Return the product of the key/values of a dictionary."""
    keys = dd.keys()
    for element in itertools.product(*dd.values()):
        yield dict(zip(keys, element))

def prep_pdbs(proteinpaths, proteinalignedpaths, boundingBoxSizes, box_padding):
     for prot in proteins:
          proteinpaths[prot] = get_protein_path(prot+".pdb")
          proteinalignedpaths[prot] = get_protein_path(prot+"_aligned"+".pdb")
          align_protein_to_inertial_axes(proteinpaths[prot],proteinalignedpaths[prot])
          ([minX, minY, minZ],[maxX, maxY, maxZ]) = get_protein_dimensions(proteinalignedpaths[prot])
          print("The Inertia Axis Aligned Bounding Box (IABB) dimensions of %s are (%.2f, %.2f, %.2f)" % (proteinalignedpaths[prot],maxX-minX, maxY-minY, maxZ-minZ))
          print("The Inertia Axis Aligned Bounding Box (IABB) volume of %s is %.2f A3" % (proteinalignedpaths[prot],(maxX-minX)*(maxY-minY)*(maxZ-minZ)))
          boundingBoxSizes[prot] = np.array([maxX-minX, maxY-minY, maxZ-minZ])
          liq_box_lengths[prot] = [u.unyt_array(boundingBoxSizes[prot]+2*box_padding, u.angstrom)]


waterModel = ["tip3"]

moleculeNameAsXML = True

# None is used as the statepoint to equilibrate the solvent
salt_strengths = [0.300, None]
#salt_strengths = [0.000, 0.075, 0.150, 0.225, 0.300, None]
proteins = [
    "6g6k",
]
sim_types = [
    "dimer", "myc_mono", "max_mono"
]
proteinpaths = dict()
proteinalignedpaths = dict()
boundingBoxSizes = dict()
liq_box_lengths = dict()

box_padding = 15
empty_space = 2
prep_pdbs(proteinpaths, proteinalignedpaths, boundingBoxSizes, box_padding+empty_space)

cations = [["SOD", "1"]]
anions = [["CLA", "-1"]]


#replicas = range(16)
#replicas = range(3)
replicas = range(3)

md_engines = ["gromacs",]
#md_engines = ["namd", "gromacs", "hoomd", "lammps-VU", "lammps-UD"]
#mc_engines = ["cassandra", "mcccs", "gomc"]
mc_engines = []
forcefields = ["a99SBdisp"]
r_cuts = [10] * u.angstrom
cutoff_styles = ["hard"]
long_range_correction = ["energy_pressure"]

masses = [18.0153] * u.amu


# 300  -> 280.0, 300.0, 320.0
temperatures = [298.15] * u.K
# 101.325, 101.325, 101.325
pressures = [101.325] * u.kPa

pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)


# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()
for prot in proteins:
    for (
        engine,
        (temp, press),
        liq_box_L,
        mass,
        rcut,
        lrc,
        cutoff_style,
        replica,
        conc,
        cat,
        an,
        wm,
        ff,
        sim_type
    ) in itertools.product(
        md_engines,
        zip(temperatures, pressures),
        liq_box_lengths[prot],
        masses,
        r_cuts,
        long_range_correction,
        cutoff_styles,
        replicas,
        salt_strengths,
        cations,
        anions,
        waterModel,
        forcefields,
        sim_types
    ):
        statepoint = {
            "waterModel": wm,
            "salt_conc": conc,
            "cat_name": cat[0],
            "cat_val": cat[1],
            "an_name": an[0],
            "an_val": an[1],
            "pdbid" : prot,
            "engine": engine,
            "replica": replica,
            "temperature": np.round(
                temp.to_value("K"),
                decimals=3,
            ).item(),
            "pressure": np.round(press.to_value("kPa"), decimals=3).item(),
            "box_L_liq_x": np.round(
                liq_box_L[0].to_value("nm"),
                decimals=3,
            ).item()
            if liq_box_L[0]
            else None,
            "box_L_liq_y": np.round(
                liq_box_L[1].to_value("nm"),
                decimals=3,
            ).item()
            if liq_box_L[1]
            else None,
            "box_L_liq_z": np.round(
                liq_box_L[2].to_value("nm"),
                decimals=3,
            ).item()
            if liq_box_L[2]
            else None,
            "mass": np.round(
                mass.to_value("amu"),
                decimals=3,
            ).item(),
            "forcefield_name": ff,
            "cutoff_style": cutoff_style,
            "long_range_correction": lrc,
            "rcut": np.round(
                rcut.to_value("nm"),
                decimals=3,
            ).item(),
            "sim_type" : sim_type
        }
        total_statepoints.append(statepoint)

for sp in total_statepoints:
    pr.open_job(
        statepoint=sp,
    ).init()
