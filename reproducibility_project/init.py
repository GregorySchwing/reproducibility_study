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

solvent_molecules = [
    "tip3p",
    "spce",
    "opc3",
    "tip4p_ew",
    "tip4p_2005",
    "tip4p_d",
    "a99SB_disp",
    "opc",
]

proteins = [
    "6g6k",
]
proteinpaths = dict()
proteinalignedpaths = dict()
boundingBoxSizes = dict()
liq_box_lengths = dict()
box_padding = 10
empty_space = 2
prep_pdbs(proteinpaths, proteinalignedpaths, boundingBoxSizes, box_padding+empty_space)
moleculeNameAsXML = True

replicas = range(1)
simulation_engines = [
    "gomc",
]
md_engines = ["gromacs", "hoomd", "lammps-VU", "lammps-UD", "namd"]
mc_engines = ["cassandra", "mcccs", "gomc"]
forcefields = {}
r_cuts = {}
cutoff_styles = ["hard"]
long_range_correction = ["energy_pressure"]

for key in solvent_molecules:
    if "UA" in key:
        if "benz" not in key:
            forcefields[key] = "trappe-ua"
        else:
            forcefields[key] = "benzene-ua"
        r_cuts[key] = 14.0 * u.angstrom
    elif "SPCE" in key:
        forcefields[key] = "spce"
        r_cuts[key] = 9 * u.angstrom
    elif moleculeNameAsXML:
        forcefields[key] = "custom"
        r_cuts[key] = 9 * u.angstrom
    else :
        forcefields[key] = "oplsaa"
        r_cuts[key] = 10 * u.angstrom
g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
masses = {
    "tip3p": [18.0153] * u.amu,
    "spce": [18.0153] * u.amu,
    "opc3": [18.0153] * u.amu,
    "tip4p_ew": [18.0153] * u.amu,
    "tip4p_2005": [18.0153] * u.amu,
    "tip4p_d": [18.0153] * u.amu,
    "a99SB_disp": [18.0153] * u.amu,
    "opc": [18.0153] * u.amu,
}
init_density_liq = {
    "tip3p": [98.0] * g_per_cm3,
    "spce": [99.3] * g_per_cm3,
    "opc3": [99.1] * g_per_cm3,
    "tip4p_ew": [99.6] * g_per_cm3,
    "tip4p_2005": [99.7] * g_per_cm3,
    "tip4p_d": [99.3] * g_per_cm3,
    "a99SB_disp": [99.6] * g_per_cm3,
    "opc": [99.7] * g_per_cm3,
}
init_density_vap = {
    "tip3p": [None],
    "spce": [None],
    "opc3": [None],
    "tip4p_ew": [None],
    "tip4p_2005": [None],
    "tip4p_d": [None],
    "a99SB_disp": [None],
    "opc": [None],
}
# 300  -> 280.0, 300.0, 320.0
temperatures = {
    "tip3p": [300.0] * u.K,
    "spce": [300.0] * u.K,
    "opc3": [300.0] * u.K,
    "tip4p_ew": [300.0] * u.K,
    "tip4p_2005": [300.0] * u.K,
    "tip4p_d": [300.0] * u.K,
    "a99SB_disp": [300.0] * u.K,
    "opc": [300.0] * u.K,
}
# 101.325, 101.325, 101.325
pressures = {
    "tip3p": [101.325] * u.kPa,
    "spce": [101.325] * u.kPa,
    "opc3": [101.325] * u.kPa,
    "tip4p_ew": [101.325] * u.kPa,
    "tip4p_2005": [101.325] * u.kPa,
    "tip4p_d": [101.325] * u.kPa,
    "a99SB_disp": [101.325] * u.kPa,
    "opc": [101.325] * u.kPa,
}
#1100, 1100, 1100
N_liq_solvent_molecules = {
    "tip3p": [1100],
    "spce": [1100],
    "opc3": [1100],
    "tip4p_ew": [1100],
    "tip4p_2005": [1100],
    "tip4p_d": [1100],
    "a99SB_disp": [1100],
    "opc": [1100],
}

N_vap_solvent_molecules = {
    "spce": [None],
    "tip3p": [None],
    "opc3": [None],
    "tip4p_ew": [None],
    "tip4p_2005": [None],
    "tip4p_d": [None],
    "a99SB_disp": [None],
    "opc": [None],
}

vap_box_lengths = {
    "tip3p": [None],
    "spce": [None],
    "opc3": [None],
    "tip4p_ew": [None],
    "tip4p_2005": [None],
    "tip4p_d": [None],
    "a99SB_disp": [None],
    "opc": [None],
}

ensembles = {
    "tip3p": ["NPT", None],
    "spce": ["NPT", None],
    "opc3": ["NPT", None],
    "tip4p_ew": ["NPT", None],
    "tip4p_2005": ["NPT", None],
    "tip4p_d": ["NPT", None],
    "a99SB_disp": ["NPT", None],
    "opc": ["NPT", None],
}


pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)
# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()
for prot in proteins:
    for molecule in solvent_molecules:
        for (
            engine,
            ensemble,
            (temp, press),
            n_liq,
            liq_box_L,
            n_vap,
            vap_box_L,
            (init_liq_den, init_vap_den),
            mass,
            lrc,
            cutoff_style,
            replica,
        ) in itertools.product(
            simulation_engines,
            ensembles[molecule],
            zip(temperatures[molecule], pressures[molecule]),
            N_liq_solvent_molecules[molecule],
            liq_box_lengths[prot],
            N_vap_solvent_molecules[molecule],
            vap_box_lengths[molecule],
            zip(init_density_liq[molecule], init_density_vap[molecule]),
            masses[molecule],
            long_range_correction,
            cutoff_styles,
            replicas,
        ):
            statepoint = {
                "molecule": molecule,
                "pdbid" : prot,
                "engine": engine,
                "replica": replica,
                "temperature": np.round(
                    temp.to_value("K"),
                    decimals=3,
                ).item(),
                "pressure": np.round(press.to_value("kPa"), decimals=3).item(),
                "ensemble": ensemble if ensemble else None,
                "N_liquid": n_liq,
                "N_vap": n_vap if n_vap else None,
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
                "box_L_vap": np.round(
                    vap_box_L.to_value("nm"),
                    decimals=3,
                ).item()
                if vap_box_L
                else None,
                "init_liq_den": np.round(
                    init_liq_den.to_value(g_per_cm3),
                    decimals=3, 
                ).item(),
                "init_vap_den": np.round(
                    init_vap_den.to_value(g_per_cm3),
                    decimals=3,
                ).item()
                if init_vap_den
                else None,
                "mass": np.round(
                    mass.to_value("amu"),
                    decimals=3,
                ).item(),
                "forcefield_name": forcefields[molecule],
                "cutoff_style": cutoff_style,
                "long_range_correction": lrc,
                "r_cut": np.round(
                    r_cuts[molecule].to_value("nm"),
                    decimals=3,
                ).item(),
            }
            total_statepoints.append(statepoint)

# print(len(total_statepoints))
indices_to_remove = set()
for i, sp in enumerate(total_statepoints):
    # filter gemc ensembles from md engines
    if sp["ensemble"] == "GEMC-NVT" and sp["engine"] in md_engines:
        indices_to_remove.add(i)

    if sp["ensemble"] == "NPT":
        sp["N_vap"] = None
        sp["box_L_vap"] = None
        sp["init_vap_den"] = None

    if sp["ensemble"] is None:
        indices_to_remove.add(i)

    if (
        sp["engine"] in mc_engines
        and sp["molecule"] == "pentaneUA-flexible_bonds"
    ):
        indices_to_remove.add(i)
    if (
        "lammps" in sp["engine"]
        and sp["molecule"] == "pentaneUA-constrain_bonds"
    ):
        indices_to_remove.add(i)


# now reverse sort the set and remove from inital list
# must be reverse sorted to remove indices on the list in place
# otherwise the list will change size and the indices would change
# print(len(indices_to_remove))
sorted_indicies_to_delete = sorted(list(indices_to_remove), reverse=True)
for idx in sorted_indicies_to_delete:
    del total_statepoints[idx]

for sp in total_statepoints:
    pr.open_job(
        statepoint=sp,
    ).init()
