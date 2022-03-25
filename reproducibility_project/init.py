"""Initialize signac statepoints."""
import itertools
import os

import numpy as np
import signac
import unyt as u
from numpy import ModuleDeprecationWarning


def dict_product(dd):
    """Return the product of the key/values of a dictionary."""
    keys = dd.keys()
    for element in itertools.product(*dd.values()):
        yield dict(zip(keys, element))


molecules = [
    "waterTIP3P",
    "waterSPCE",
    "waterTI4PEW",
    "waterTI4P2005",
    "waterTI4PD",
    "waterA99SBdisp",
    "waterOPC",
]
replicas = range(16)
simulation_engines = [
    "cassandra",
    "mcccs",
    "gomc",
    "gromacs",
    "hoomd",
    "lammps-VU",
    "lammps-UD",
]
md_engines = ["gromacs", "hoomd", "lammps-VU", "lammps-UD"]
mc_engines = ["cassandra", "mcccs", "gomc"]
forcefields = {}
r_cuts = {}
cutoff_styles = ["hard"]
long_range_correction = ["energy_pressure"]
for key in molecules:
    if "UA" in key:
        if "benz" not in key:
            forcefields[key] = "trappe-ua"
        else:
            forcefields[key] = "benzene-ua"
        r_cuts[key] = 14.0 * u.angstrom
    elif "SPCE" in key:
        forcefields[key] = "spce"
        r_cuts[key] = 9 * u.angstrom
    else:
        forcefields[key] = "oplsaa"
        r_cuts[key] = 10 * u.angstrom
g_per_cm3 = u.g / (u.cm * u.cm * u.cm)
masses = {
    "waterTIP3P": [18.0153] * u.amu,
    "waterSPCE": [18.0153] * u.amu,
    "waterTI4PEW": [18.0153] * u.amu,
    "waterTI4P2005": [18.0153] * u.amu,
    "waterTI4PD": [18.0153] * u.amu,
    "waterA99SBdisp": [18.0153] * u.amu,
    "waterOPC": [18.0153] * u.amu,
}
init_density_liq = {
    "waterTIP3P": [0.98] * g_per_cm3,
    "waterSPCE": [0.993] * g_per_cm3,
    "waterTI4PEW": [0.996] * g_per_cm3,
    "waterTI4P2005": [0.997] * g_per_cm3,
    "waterTI4PD": [0.993] * g_per_cm3,
    "waterA99SBdisp": [0.996] * g_per_cm3,
    "waterOPC": [0.997] * g_per_cm3,
}
init_density_vap = {
    "waterTIP3P": [None],
    "waterSPCE": [None],
    "waterTI4PEW": [None],
    "waterTI4P2005": [None],
    "waterTI4PD": [None],
    "waterA99SBdisp": [None],
    "waterOPC": [None],
}
temperatures = {
    "waterTIP3P": [280.0, 300.0, 320.0] * u.K,
    "waterSPCE": [280.0, 300.0, 320.0] * u.K,
    "waterTI4PEW": [280.0, 300.0, 320.0] * u.K,
    "waterTI4P2005": [280.0, 300.0, 320.0] * u.K,
    "waterTI4PD": [280.0, 300.0, 320.0] * u.K,
    "waterA99SBdisp": [280.0, 300.0, 320.0] * u.K,
    "waterOPC": [280.0, 300.0, 320.0] * u.K,
}

pressures = {
    "waterTIP3P": [101.325, 101.325, 101.325] * u.kPa,
    "waterSPCE": [101.325, 101.325, 101.325] * u.kPa,
    "waterTI4PEW": [101.325, 101.325, 101.325] * u.kPa,
    "waterTI4P2005": [101.325, 101.325, 101.325] * u.kPa,
    "waterTI4PD": [101.325, 101.325, 101.325] * u.kPa,
    "waterA99SBdisp": [101.325, 101.325, 101.325] * u.kPa,
    "waterOPC": [101.325, 101.325, 101.325] * u.kPa,
}

N_liq_molecules = {
    "waterSPCE": [1100, 1100, 1100],
    "waterTIP3P": [1100, 1100, 1100],
    "waterSPCE": [1100, 1100, 1100],
    "waterTI4PEW": [1100, 1100, 1100],
    "waterTI4P2005": [1100, 1100, 1100],
    "waterTI4PD": [1100, 1100, 1100],
    "waterA99SBdisp": [1100, 1100, 1100],
    "waterOPC": [1100, 1100, 1100],
}

N_vap_molecules = {
    "waterSPCE": [None],
    "waterTIP3P": [None],
    "waterSPCE": [None],
    "waterTI4PEW": [None],
    "waterTI4P2005": [None],
    "waterTI4PD": [None],
    "waterA99SBdisp": [None],
    "waterOPC": [None],
}

liq_box_lengths = {
    "waterTIP3P": [32.07] * u.angstrom,
    "waterSPCE": [32.07] * u.angstrom,
    "waterTI4PEW": [32.07] * u.angstrom,
    "waterTI4P2005": [32.07] * u.angstrom,
    "waterTI4PD": [32.07] * u.angstrom,
    "waterA99SBdisp": [32.07] * u.angstrom,
    "waterOPC": [32.07] * u.angstrom,
}
}

vap_box_lengths = {
    "waterSPCE": [None],
    "waterTIP3P": [None],
    "waterSPCE": [None],
    "waterTI4PEW": [None],
    "waterTI4P2005": [None],
    "waterTI4PD": [None],
    "waterA99SBdisp": [None],
    "waterOPC": [None],
}

ensembles = {
    "waterSPCE": ["NPT", None],
    "waterSPCE": ["NPT", None],
    "waterTIP3P": ["NPT", None],
    "waterSPCE": ["NPT", None],
    "waterTI4PEW": ["NPT", None],
    "waterTI4P2005": ["NPT", None],
    "waterTI4PD": ["NPT", None],
    "waterA99SBdisp": ["NPT", None],
    "waterOPC": ["NPT", None],
}


pr_root = os.path.join(os.getcwd(), "src")
pr = signac.get_project(pr_root)

# ignore statepoints that are not being tested (gemc only for methane, pentane)
# filter the list of dictionaries
total_statepoints = list()
for molecule in molecules:
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
        N_liq_molecules[molecule],
        liq_box_lengths[molecule],
        N_vap_molecules[molecule],
        vap_box_lengths[molecule],
        zip(init_density_liq[molecule], init_density_vap[molecule]),
        masses[molecule],
        long_range_correction,
        cutoff_styles,
        replicas,
    ):
        statepoint = {
            "molecule": molecule,
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
            "box_L_liq": np.round(
                liq_box_L.to_value("nm"),
                decimals=3,
            ).item()
            if liq_box_L
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
