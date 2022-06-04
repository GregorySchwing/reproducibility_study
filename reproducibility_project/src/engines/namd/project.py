"""GOMC's setup for signac, signac-flow, signac-dashboard for this study."""
# project.py

import os
import subprocess

import flow
import matplotlib.pyplot as plt
import math
# from flow.environment import StandardEnvironment
import mbuild as mb
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control
import numpy as np
import pandas as pd
import pymbar
import signac
import unyt as u
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment
from flow.environments.incite import SummitEnvironment
from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.molecules.system_builder import (
    construct_system,
)
from src.proteins.protein_processor import (
    align_protein_to_inertial_axes,
    get_protein_dimensions,
    get_protein_path,
    merge_solv_and_solute,
    fix_segment,
    ionize
)

from reproducibility_project.src.utils.forcefields import get_ff_path
from reproducibility_project.src.utils.plotting import plot_data_with_t0_line


class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()

class Grid(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
    #Subclass of DefaultSlurmEnvironment for WSU's Grid cluster.

    hostname_pattern = r".*\.grid\.wayne\.edu"
    template = "grid.sh"


#equilibrateSolvent = Project.make_group(name="equilibrateSolvent")
#prepareProteinSimulation = Project.make_group(name="prepareProteinSimulation")

# ******************************************************
# users typical variables, but not all (start)
# ******************************************************
# set binary path to gomc binary files (the bin folder).
# If the gomc binary files are callable directly from the terminal without a path,
# please just enter and empty string (i.e., "" or '')
# namd_binary_path = "/wsu/home/hf/hf68/hf6839/GOMC_dev_2_22_22/bin"

namd_binary_path = os.path.join(os.getcwd(), "bin")
namd2_binary = "namd2"
namd3_binary = "namd3"


# number of MC cycles
MC_cycles_melt_equilb_NVT = 5 * 10 ** 3  # set value for paper = 5 * 10 ** 3
MC_cycles_equilb_NVT = 5 * 10 ** 3  # set value for paper = 5 * 10 ** 3
MC_cycles_equilb_design_ensemble = (
    40 * 10 ** 3
)  # set value for paper = 40 * 10 ** 3
MC_cycles_production = 120 * 10 ** 3  # set value for paper = 120 * 10 ** 3

output_data_every_X_MC_cycles = 10  # set value for paper = 10

# max number of equilibrium selected runs
equilb_design_ensemble_max_number = 3

# force field (FF) file for all simulations in that job
# Note: do not add extensions
ff_filename_str = "in_FF"
"""
min_steps = 500
nvt_eq_steps = 5000
npt_eq_steps = 5000
production_steps = 10000
single_production_run_steps = 5000
"""
min_steps = 10000
nvt_eq_steps = 250000
npt_eq_steps = 1000000
production_steps = 5000000
single_production_run_steps = 500000

num_cycles = math.ceil(production_steps / single_production_run_steps)
# initial mosdef structure and coordinates
# Note: do not add extensions
mosdef_structure_box_0_name_str = "mosdef_box_0"
mosdef_structure_box_1_name_str = "mosdef_box_1"

# melt equilb simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
min_NVT_control_file_name_str = "em"

# equilb simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
equilb_NVT_control_file_name_str = "nvt_eq"

# equilb simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
equilb_NPT_control_file_name_str = "npt_eq"

# The production run using the ensemble used for the simulation design, which
# includes the simulation runs GOMC control file input and simulation outputs
# Note: do not add extensions
production_control_file_name_str = "npt_prod"


path_from_job_to_box_inputs = "../../"

walltime_mosdef_hr = 4
walltime_gomc_hr = 4
memory_needed = 16

use_pymbar = True  # True of False

ff_info_dict = {
    "trappe-ua": {
        "ngpu": 0,
        "ncpu": 1,
        "Ewald": False,
        "ElectroStatic": False,
        "VDWGeometricSigma": False,
    },
    "benzene-ua": {
        "ngpu": 0,
        "ncpu": 1,
        "Ewald": False,
        "ElectroStatic": False,
        "VDWGeometricSigma": False,
    },
    "spce": {
        "ngpu": 1,
        "ncpu": 4,
        "Ewald": True,
        "ElectroStatic": True,
        "VDWGeometricSigma": False,
    },
    "oplsaa": {
        "ngpu": 0,
        "ncpu": 8,
        "Ewald": True,
        "ElectroStatic": True,
        "VDWGeometricSigma": True,
    },
    "custom": {
        "ngpu": 1,
        "ncpu": 4,
        "Ewald": True,
        "ElectroStatic": True,
        "VDWGeometricSigma": False,
    },
}

# ******************************************************
# users typical variables, but not all (end)
# ******************************************************


# ******************************************************
# signac and GOMC-MOSDEF code (start)
# ******************************************************
# checking if the GOMC control file is written for the melt equilb NVT run

# ******************************************************
# ******************************************************
# create some initial variable to be store in each jobs
# directory in an additional json file, and test
# to see if they are written (start).
# ******************************************************
# ******************************************************
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
def part_1a_initial_data_input_to_json(job):
    """Check that the initial job data is written to the json files."""
    data_written_bool = False
    if job.isfile(f"{'signac_job_document.json'}"):
        data_written_bool = True

    return data_written_bool

#@equilibrateSolvent
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(lambda j: j.sp.replica == 0)
@Project.pre(lambda j: j.sp.salt_conc == None)
@Project.post(part_1a_initial_data_input_to_json)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get("ncpu"),
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def initial_parameters(job):
    """Set the initial job parameters into the jobs doc json file."""
    job.doc.cycle = 0
    job.doc.num_cycles = num_cycles
    # select
    job.doc.ngpu = ff_info_dict.get(job.sp.forcefield_name).get("ngpu")
    if job.doc.ngpu == 0:
        job.doc.cpu_or_gpu = "CPU"
    elif job.doc.ngpu == 1:
        job.doc.cpu_or_gpu = "GPU"
    else:
        raise ValueError(
            "CPU and GPU can not be deterimined as force field (FF) is not available in the selection."
        )

    # FF type to directory and path
    job.doc.forcefield_directory_name = get_ff_path(job.sp.forcefield_name, job.sp.molecule)

    # reformat ensembles for input to GOMC
    if job.sp.ensemble in ["GEMC-NPT", "GEMC_NPT"]:
        job.doc.production_ensemble = "GEMC_NPT"
    elif job.sp.ensemble in ["GEMC-NVT", "GEMC_NVT"]:
        job.doc.production_ensemble = "GEMC_NVT"
    else:
        job.doc.production_ensemble = job.sp.ensemble

    # list replica seed numbers
    replica_no_to_seed_dict = {
        0: 0,
        1: 1,
        2: 2,
        3: 3,
        4: 4,
        5: 5,
        6: 6,
        7: 7,
        8: 8,
        9: 9,
        10: 10,
        11: 11,
        12: 12,
        13: 13,
        14: 14,
        15: 15,
        16: 16,
        17: 17,
        18: 18,
        19: 19,
        20: 20,
    }

    job.doc.replica_number_int = replica_no_to_seed_dict.get(
        int(job.sp.replica)
    )

    # set rcut, ewalds
    job.doc.Ewald = ff_info_dict.get(job.sp.forcefield_name).get("Ewald")
    job.doc.ElectroStatic = ff_info_dict.get(job.sp.forcefield_name).get(
        "ElectroStatic"
    )
    job.doc.VDWGeometricSigma = ff_info_dict.get(job.sp.forcefield_name).get(
        "VDWGeometricSigma"
    )

    # set the initial iteration number of the simulation
    job.doc.equilb_design_ensemble_dict = {}
    job.doc.equilb_design_ensemble_number = 0
    job.doc.equilb_design_ensemble_max_number = (
        equilb_design_ensemble_max_number
    )
    job.doc.equilb_design_ensemble_max_number_under_limit = True
    job.doc.stable_equilb_design_ensemble = False

    job.doc.production_run_ensemble_dict = {}

    if job.doc.production_ensemble in ["NVT", "NPT"]:
        job.doc.N_liquid = job.sp.N_liquid
        job.doc.N_vap = 0
        job.doc.N_total = job.sp.N_liquid
    elif job.doc.production_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
        job.doc.N_liquid = job.sp.N_liquid
        job.doc.N_vap = job.sp.N_vap
        job.doc.N_total = job.sp.N_liquid + job.sp.N_vap
    # reject if set design ensemble is NVT
    # currently pymbar is done of of density, which will not work for NVT
    if job.doc.production_ensemble == "NVT":
        raise ValueError(
            "ERROR: The NVT ensemble is not currently available for this project.py "
            "script, as pymbar is done based on density, which will not work for NVT"
        )

    # select binary path and binary file
    job.doc.namd_binary_path = namd_binary_path

    if job.doc.production_ensemble in ["NPT", "NVT"]:
        job.doc.melt_NVT_gomc_binary_file = namd2_binary
        job.doc.equilb_NVT_gomc_binary_file = namd3_binary
        job.doc.equilb_NPT_gomc_binary_file = namd3_binary
        job.doc.equilb_design_ensemble_gomc_binary_file = (
            namd2_binary
        )
    else:
        raise ValueError(
            "ERROR: A wrong ensemble has been specified for the namd binary file"
        )

    if job.doc.production_ensemble in ["NPT"]:
        job.doc.production_ensemble_gomc_binary_file = (
            namd3_binary
        )
    elif job.doc.production_ensemble in ["NVT"]:
        job.doc.production_ensemble_gomc_binary_file = (
            namd3_binary
        )
    else:
        raise ValueError(
            "ERROR: A wrong ensemble has been specified for the gomc binary file"
        )


# ******************************************************
# ******************************************************
# create some initial variable to be store in each jobs
# directory in an additional json file, and test
# to see if they are written (end).
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC psf, pdb, and force field (FF) files were written (start)
# ******************************************************
# ******************************************************


@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
def part_1b_under_equilb_design_ensemble_run_limit(job):
    """Check that the equilbrium design ensemble run is under it's run limit."""
    try:
        if (
            job.doc.equilb_design_ensemble_number
            >= job.doc.equilb_design_ensemble_max_number
        ):
            job.doc.equilb_design_ensemble_max_number_under_limit = False
            return job.doc.equilb_design_ensemble_max_number_under_limit

        else:
            return True
    except:
        return False




# ******************************************************
# ******************************************************
# check if GOMC psf, pdb, and FF files were written (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC control file was written (start)
# ******************************************************
# ******************************************************
# function for checking if the GOMC control file is written
def gomc_control_file_written(job, control_filename_str):
    """General check that the gomc control files are written."""
    file_written_bool = False
    control_file = f"{control_filename_str}.conf"

    return job.isfile(control_file)


# checking if the GOMC control file is written for the melt equilb NVT run
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_2a_min_NVT_control_file_written(job):
    """General check that the min control file is written."""
    return gomc_control_file_written(job, min_NVT_control_file_name_str)


# checking if the GOMC control file is written for the equilb NVT run
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_2b_equilb_NVT_control_file_written(job):
    """General check that the equilb_NVT_control (run temperature) gomc control file is written."""
    return gomc_control_file_written(job, equilb_NVT_control_file_name_str)


# checking if the GOMC control file is written for the equilb run with the selected ensemble
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_2c_equilb_NPT_control_file_written(job):
    """General check that the equilb_NPT_control (run temperature) gomc control file is written."""
    return gomc_control_file_written(job, equilb_NPT_control_file_name_str)


# checking if the GOMC control file is written for the production run
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(part_1a_initial_data_input_to_json)
@flow.with_job
def part_2d_production_control_file_written(job):
    """General check that the prod_NPT_control (run temperature) gomc control file is written."""
    try:
        cycleList = list(range(0, job.doc.num_cycles))
        allConfsExist = True
        for cycle in cycleList:
            cycleExists = gomc_control_file_written(job, production_control_file_name_str+"_"+str(job.doc.cycle))
            allConfsExist = allConfsExist and cycleExists
        return allConfsExist
    except:
        return False

# check if GOMC-MOSDEF wrote the gomc files
# @Project.pre(select_production_ensemble)
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def mosdef_input_written(job):
    """Check that the initial job data is written to the json files."""
    data_written_bool = False
    if job.isfile(f"{'signac_job_document.json'}"):
        data_written_bool = True

    return data_written_bool


    """Check that the mosdef files (psf, pdb, and force field (FF) files) are written .
    file_written_bool = False
    try:
        if job.sp.ensemble in ["NPT", "NVT"]:
            if (
                job.isfile(f"{path_from_job_to_box_inputs}/{ff_filename_str}.inp")
                and job.isfile(
                    f"{path_from_job_to_box_inputs}/{mosdef_structure_box_0_name_str}.psf"
                )
                and job.isfile(
                    f"{path_from_job_to_box_inputs}/{mosdef_structure_box_0_name_str}.pdb"
                )
            ):
                file_written_bool = True
        elif job.sp.ensemble in ["GCMC", "GEMC_NPT", "GEMC_NPT"]:
            if (
                job.isfile(f"{path_from_job_to_box_inputs}/{ff_filename_str}.inp")
                and job.isfile(
                    f"{path_from_job_to_box_inputs}/{mosdef_structure_box_0_name_str}.psf"
                )
                and job.isfile(
                    f"{path_from_job_to_box_inputs}/{mosdef_structure_box_0_name_str}.pdb"
                )
                and job.isfile(
                    f"{path_from_job_to_box_inputs}/{mosdef_structure_box_1_name_str}.psf"
                )
                and job.isfile(
                    f"{path_from_job_to_box_inputs}/{mosdef_structure_box_1_name_str}.pdb"
                )
            ):
                file_written_bool = True
        return file_written_bool
    except:
        return False
    """

@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_2a_solvated(job):
    """Check that the initial job data is written to the json files."""
    data_written_bool = False
    saltless_sp = job.statepoint()
    saltless_sp['salt_conc']=None
    saltless_sp['replica']=0
    #print("statepoint desalted",saltless_sp)
    #print(Project.doc)
    print(Project.find_jobs())
    #res = Project.find_jobs(saltless_sp)
    #jobs = list(Project.find_jobs(saltless_sp))
    #for job in jobs:
    #    print(job.fn("solvated.pdb"))
    #print(res)
    #for job in res:
    #    print("job id", Project.open_job(job).id)
    #if res.next().isfile(f"{'solvated.pdb'}"):
    #    data_written_bool = True

    return data_written_bool


@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_2a_ionized(job):
    """Check that the initial job data is written to the json files."""
    data_written_bool = False
    if job.isfile(f"{'ionized.pdb'}"):
        data_written_bool = True

    return data_written_bool

# ******************************************************
# ******************************************************
# check if GOMC control file was written (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC simulations started (start)
# ******************************************************
# ******************************************************
# function for checking if GOMC simulations are started
def gomc_simulation_started(job, control_filename_str):
    """General check to see if the gomc simulation is started."""
    output_started_bool = False
    if job.isfile("out_{}.dat".format(control_filename_str)):
        output_started_bool = True

    return output_started_bool


# check if melt equilb_NVT GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(lambda j: j.sp.replica == 0)
@flow.with_job
def part_3a_output_melt_equilb_NVT_started(job):
    """Check to see if the melt_equilb_NVT (high temperature) gomc simulation is started."""
    return gomc_simulation_started(job, min_NVT_control_file_name_str)


# check if equilb_NVT GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_3b_output_equilb_NVT_started(job):
    """Check to see if the equilb_NVT (set temperature) gomc simulation is started."""
    return gomc_simulation_started(job, equilb_NVT_control_file_name_str)

# check if equilb_NVT GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_3c_output_equilb_NPT_started(job):
    """Check to see if the equilb_NVT (set temperature) gomc simulation is started."""
    return gomc_simulation_started(job, equilb_NPT_control_file_name_str)

# check if production GOMC run is started by seeing if the GOMC consol file and the merged psf exist
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(part_1a_initial_data_input_to_json)
@flow.with_job
def part_3d_output_production_run_started(job):
    """Check to see if the production run (set temperature) gomc simulation is started."""
    try:
        return gomc_simulation_started(job, production_control_file_name_str+"_"+str(job.doc.cycle))

    except:
        return False


# ******************************************************
# ******************************************************
# check if GOMC simulations started (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# check if GOMC simulation are completed properly (start)
# ******************************************************
# ******************************************************

# function for checking if GOMC simulations are completed properly
def gomc_sim_completed_properly(job, control_filename_str):
    """General check to see if the gomc simulation was completed properly."""
    job_run_properly_bool = False
    output_log_file = "out_{}.dat".format(control_filename_str)

    if job.isfile(output_log_file):
        with open(f"{output_log_file}", 'rb') as f:
            try:  # catch OSError in case of a one line file 
                f.seek(-20, os.SEEK_END)
                while f.read(1) != b'\n':
                    f.seek(-2, os.SEEK_CUR)
            except OSError:
                f.seek(0)
            last_line = f.readline().decode()
            if("End of program" in last_line):
                job_run_properly_bool = True
            else:
                job_run_properly_bool = False

    return job_run_properly_bool


# check if melt equilb NVT GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_4a_job_min_NVT_completed_properly(job):
    """Check to see if the melt_equilb_NVT (high temperature) gomc simulation was completed properly."""
    return gomc_sim_completed_properly(
        job, min_NVT_control_file_name_str
    )


# check if equilb NVT GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_4b_job_equilb_NVT_completed_properly(job):
    """Check to see if the equilb_NVT (set temperature) gomc simulation was completed properly."""
    return gomc_sim_completed_properly(job, equilb_NVT_control_file_name_str)

@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@flow.with_job
def part_4c_job_equilb_NPT_completed_properly(job):
    """Check to see if the equilb_NVT (set temperature) gomc simulation was completed properly."""
    return gomc_sim_completed_properly(job, equilb_NPT_control_file_name_str)


# check if production GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(part_1a_initial_data_input_to_json)
@flow.with_job
def part_4d_job_production_run_completed_properly(job):
    """Check to see if the production run (set temperature) gomc simulation was completed properly."""
    try:
        return gomc_sim_completed_properly(job, production_control_file_name_str+"_"+str(job.doc.cycle))
    except:
        return False

# check if production GOMC run completed by checking the end of the GOMC consol file
@Project.label
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(part_1a_initial_data_input_to_json)
@flow.with_job
def part_4e_job_final_production_run_completed_properly(job):
    """Check to see if the production run (set temperature) gomc simulation was completed properly."""
    try:
        return job.doc.cycle == job.doc.num_cycles
    except:
        return False

# ******************************************************
# ******************************************************
# check if GOMC simulation are completed properly (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# build system, with option to write the force field (force field (FF)), pdb, psf files.
# Note: this is needed to write GOMC control file, even if a restart (start)
# ******************************************************
# build system
def build_charmm(job, write_files=True):
    """Build the Charmm object and potentially write the pdb, psd, and force field (FF) files."""
    print("#**********************")
    print("Started: GOMC Charmm Object")
    print("#**********************")

    print("Running: box packing")
    [box_0, box_1] = construct_system(job.sp)
    print("Completed: box packing")

    Bead_to_atom_name_dict = {
        "_CH4": "C",
        "_CH3": "C",
        "_CH2": "C",
        "_CH": "C",
        "_HC": "C",
    }

    # Name all water models WAT for vmd recongnition these are waters.
    #Molecule_ResName_List = [job.sp.molecule]
    Molecule_ResName_List = ["WAT"]

    if job.sp.molecule in ["waterSPCE", "benzeneUA"]:
        gomc_fix_bonds_angles_list = Molecule_ResName_List
    else:
        gomc_fix_bonds_angles_list = ["WAT"]

    if job.doc.production_ensemble in ["NVT", "NPT"]:
        charmm = mf_charmm.Charmm(
            box_0,
            mosdef_structure_box_0_name_str,
            structure_box_1=None,
            filename_box_1=None,
            ff_filename=ff_filename_str,
            forcefield_selection=job.doc.forcefield_directory_name,
            residues=Molecule_ResName_List,
            bead_to_atom_name_dict=Bead_to_atom_name_dict,
            gomc_fix_bonds_angles=None,
        )

    elif job.doc.production_ensemble in ["GCMC", "GEMC_NVT", "GEMC_NPT"]:
        charmm = mf_charmm.Charmm(
            box_0,
            mosdef_structure_box_0_name_str,
            structure_box_1=box_1,
            filename_box_1=mosdef_structure_box_1_name_str,
            ff_filename=ff_filename_str,
            forcefield_selection=job.doc.forcefield_directory_name,
            residues=Molecule_ResName_List,
            bead_to_atom_name_dict=Bead_to_atom_name_dict,
            gomc_fix_bonds_angles=None,
        )

    if write_files == True:
        charmm.write_inp()

        charmm.write_psf()

        charmm.write_pdb()

    print("#**********************")
    print("Completed: GOMC Charmm Object")
    print("#**********************")

    return charmm



# ******************************************************
# ******************************************************
# build system, with option to write the force field (FF), pdb, psf files.
# Note: this is needed to write GOMC control file, even if a restart (end)
# ******************************************************

# ******************************************************
# ******************************************************
# Creating GOMC files (pdb, psf, force field (FF), and gomc control files (start)
# ******************************************************
# ******************************************************
#@equilibrateSolvent
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(lambda j: j.sp.replica == 0)
@Project.pre(part_1a_initial_data_input_to_json)
@Project.pre(part_1b_under_equilb_design_ensemble_run_limit)
@Project.post(part_2a_min_NVT_control_file_written)
@Project.post(part_2b_equilb_NVT_control_file_written)
@Project.post(part_2c_equilb_NPT_control_file_written)
@Project.post(part_2d_production_control_file_written)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def build_psf_pdb_ff_gomc_conf(job):

    from reproducibility_project.src.engine_input.namd import conf

    """Build the Charmm object and write the pdb, psf, and force field (FF) files for all the simulations in the workspace."""
    charmm_object_with_files = build_charmm(job, write_files=True)

    # Modify mdp files according to job statepoint parameters
    cutoff_styles = {"hard": "None", "shift": "Potential-shift"}
    lrcs = {"None": "no", "energy_pressure": "EnerPres"}

    pressure = job.sp.pressure * u.kPa
    conf_abs_path = os.path.dirname(os.path.abspath(conf.__file__))
    mdps = {
        "em": {
            "fname": "em.conf",
            "template": f"{conf_abs_path}/melt_template_protein.inp.jinja",
            "water-template": f"{conf_abs_path}/melt_template_water.inp.jinja",
            "data": {
                "is_prod" : False,
                "cycle" : job.doc.cycle,
                "structure" : job.fn("mosdef_box_0.psf"),
                "coordinates" : job.fn("mosdef_box_0.pdb"),
                "waterModel" : job.sp.waterModel,
                "parameters" : job.fn("in_FF.inp"),
                "outputname" : "em",
                "X_DIM_BOX" : job.sp.box_L_liq_x*10,
                "Y_DIM_BOX" : job.sp.box_L_liq_y*10,
                "Z_DIM_BOX" : job.sp.box_L_liq_z*10,
                "X_ORIGIN_BOX" : job.sp.box_L_liq_x/2*10,
                "Y_ORIGIN_BOX" : job.sp.box_L_liq_y/2*10,
                "Z_ORIGIN_BOX" : job.sp.box_L_liq_z/2*10,
                "cutoff": job.sp.r_cut*10,
                "switchdist": job.sp.r_cut*10 - 2,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "temp": job.sp.temperature,
                "replica": job.sp.replica,
                "lrc": lrcs[job.sp.long_range_correction],
                "minimize_steps" : min_steps,
            },
        },
        "nvt_eq": {
            "fname": "nvt_eq.conf",
            "template": f"{conf_abs_path}/nvt_eq_template_protein.inp.jinja",
            "water-template": f"{conf_abs_path}/nvt_eq_template_water.inp.jinja",
            "data": {
                "is_prod" : False,
                "cycle" : job.doc.cycle,
                "structure" : job.fn("mosdef_box_0.psf"),
                "coordinates" : job.fn("mosdef_box_0.pdb"),
                "waterModel" : job.sp.waterModel,
                "parameters" : job.fn("in_FF.inp"),
                "binary_coordinates" : job.fn("em.restart.coor"),
                "binary_boxsize" : job.fn("em.restart.xsc"),
                "binary_velocities" : job.fn("em.restart.vel"),
                "outputname" : "nvt_eq",
                "cutoff": job.sp.r_cut*10,
                "switchdist": job.sp.r_cut*10 - 2,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "temp": job.sp.temperature,
                "replica": job.sp.replica,
                "lrc": lrcs[job.sp.long_range_correction],
                "run_steps" : nvt_eq_steps,
            },
        },
        "npt_eq": {
            "fname": "npt_eq.conf",
            "template": f"{conf_abs_path}/npt_template_protein_prod.inp.jinja",
            "water-template": f"{conf_abs_path}/npt_template_water.inp.jinja",
            "data": {
                "is_prod" : False,
                "cycle" : job.doc.cycle,
                "structure" : job.fn("mosdef_box_0.psf"),
                "coordinates" : job.fn("mosdef_box_0.pdb"),
                "waterModel" : job.sp.waterModel,
                "parameters" : job.fn("in_FF.inp"),
                "binary_coordinates" : job.fn("nvt_eq.restart.coor"),
                "binary_boxsize" : job.fn("nvt_eq.restart.xsc"),
                "binary_velocities" : job.fn("nvt_eq.restart.vel"),
                "outputname" : "npt_eq",
                "cutoff": job.sp.r_cut*10,
                "switchdist": job.sp.r_cut*10 - 2,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "temp": job.sp.temperature,
                "replica": job.sp.replica,
                "lrc": lrcs[job.sp.long_range_correction],
                "run_steps" : npt_eq_steps,
            },
        },
        "npt_prod": {
            "fname": "npt_prod_",
            "template": f"{conf_abs_path}/npt_template.inp.jinja",
            "water-template": f"{conf_abs_path}/npt_template_water.inp.jinja",
            "data": {
                "is_prod" : True,
                "cycle" : job.doc.cycle,
                "structure" : job.fn("mosdef_box_0.psf"),
                "coordinates" : job.fn("mosdef_box_0.pdb"),
                "waterModel" : job.sp.waterModel,
                "parameters" : job.fn("in_FF.inp"),
                "binary_coordinates" : job.fn("npt_eq.restart.coor"),
                "binary_boxsize" : job.fn("npt_eq.restart.xsc"),
                "binary_velocities" : job.fn("npt_eq.restart.vel"),
                "outputname" : "npt_prod_"+str(job.doc.cycle),
                "cutoff": job.sp.r_cut*10,
                "switchdist": job.sp.r_cut*10 - 2,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "temp": job.sp.temperature,
                "replica": job.sp.replica,
                "lrc": lrcs[job.sp.long_range_correction],
                "run_steps" : single_production_run_steps,
            }
        }
    }

    for op, mdp in mdps.items():
        if job.sp.molecule == "waterSPCE":
            _setup_conf(
                job,
                fname=mdp["fname"],
                template=mdp["water-template"],
                data=mdp["data"],
                overwrite=True,
            )
        else:
            if (mdp["data"]["is_prod"]):
                cycle = 0
                num_cycles = job.doc.num_cycles
                while(cycle < num_cycles):
                    if(cycle > 0):
                        mdp["data"]["binary_coordinates"] = job.fn("npt_prod_"+str(cycle-1)+".restart.coor")
                        mdp["data"]["binary_boxsize"] = job.fn("npt_prod_"+str(cycle-1)+".restart.xsc")
                        mdp["data"]["binary_velocities"] = job.fn("npt_prod_"+str(cycle-1)+".restart.vel")

                    _setup_conf(
                        job,
                        fname=mdp["fname"] + str(job.doc.cycle)+".conf",
                        template=mdp["water-template"],
                        data=mdp["data"],
                        overwrite=True,
                    )
                    cycle = cycle + 1
                    job.doc.cycle = cycle
                    mdp["data"]["cycle"] = cycle
                    mdp["data"]["outputname"] = "npt_prod_"+str(cycle)
                job.doc.cycle = 0
            else:
                _setup_conf(
                    job,
                    fname=mdp["fname"],
                    template=mdp["water-template"],
                    data=mdp["data"],
                    overwrite=True,
                )





# ******************************************************
# ******************************************************
# Creating GOMC files (pdb, psf, force field (FF), and gomc control files (start)
# ******************************************************
# ******************************************************
#@prepareProteinSimulation
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(lambda j: j.sp.replica == 0)
@Project.pre(lambda j: j.sp.salt_conc == None)
@Project.pre(mosdef_input_written)
#@Project.pre(part_4c_job_equilb_NPT_completed_properly)
@Project.post(part_2a_solvated)
@Project.operation.with_directives(
    {
        "np": 1,
        "ngpu": 0,
        "memory": memory_needed,
        "walltime": walltime_mosdef_hr,
    }
)
@flow.with_job
def solvate_protein(job):
    print("#**********************")
    print("# Started the solvation process.")
    print("#**********************")
    if job.sp.pdbid:
        merge_solv_and_solute(job)
        fix_segment(job)
    else: None,

#@prepareProteinSimulation
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(lambda j: j.sp.replica == 0)
@Project.pre(lambda j: j.sp.salt_conc != None)
@Project.pre(part_2a_solvated)
@Project.post(part_2a_ionized)
@Project.operation
@flow.with_job
def ionize_protein(job):
    print("#**********************")
    print("# Started the ionization process.")
    print("#**********************")
    if job.sp.pdbid:
        ionize(job)
    else: None,

# ******************************************************
# ******************************************************
# Creating GOMC files (pdb, psf, force field (FF), and gomc control files (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# melt_NVT -starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
#@equilibrateSolvent
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(lambda j: j.sp.replica == 0)
@Project.pre(part_2a_ionized)
@Project.pre(part_2a_min_NVT_control_file_written)
@Project.post(part_3a_output_melt_equilb_NVT_started)
@Project.post(part_4a_job_min_NVT_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ncpu", 1
        ),
        "ngpu": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ngpu", 0
        ),
        "memory": memory_needed,
        "walltime": walltime_gomc_hr,
    }
)
@flow.with_job
@flow.cmd
def run_min_NVT_gomc_command(job):
    """Run the gomc melt_equilb_NVT simulation."""
    print("#**********************")
    print("# Started the run_min_NVT_gomc_command.")
    print("#**********************")

    control_file_name_str = min_NVT_control_file_name_str

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(job.doc.namd_binary_path),
        str(job.doc.melt_NVT_gomc_binary_file),
        str(ff_info_dict.get(job.sp.forcefield_name).get("ncpu")),
        str(control_file_name_str),
        str(control_file_name_str),
    )

    return run_command


# ******************************************************
# ******************************************************
# melt_NVT - including GOMC control file writing and starting the GOMC simulation (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# equilb_NVT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
#@equilibrateSolvent
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(lambda j: j.sp.replica == 0)
@Project.pre(part_2a_ionized)
@Project.pre(part_4a_job_min_NVT_completed_properly)
@Project.pre(part_2b_equilb_NVT_control_file_written)
@Project.post(part_3b_output_equilb_NVT_started)
@Project.post(part_4b_job_equilb_NVT_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get("ncpu"),
        "ngpu": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ngpu", 0
        ),
        "memory": memory_needed,
        "walltime": walltime_gomc_hr,
    }
)
@flow.with_job
@flow.cmd
def run_equilb_NVT_gomc_command(job):
    """Run the gomc equilb_NVT simulation."""
    print("#**********************")
    print("# Started the run_equilb_NVT_gomc_command.")
    print("#**********************")

    control_file_name_str = equilb_NVT_control_file_name_str

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(job.doc.namd_binary_path),
        str(job.doc.equilb_NVT_gomc_binary_file),
        str(ff_info_dict.get(job.sp.forcefield_name).get("ncpu")),
        str(control_file_name_str),
        str(control_file_name_str),
    )

    return run_command

# ******************************************************
# ******************************************************
# equilb_NPT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
# ******************************************************
# ******************************************************
# equilb NPT or GEMC-NVT - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
#@equilibrateSolvent
@Project.pre(lambda j: j.sp.engine == "namd")
@Project.pre(lambda j: j.sp.replica == 0)
@Project.pre(part_2a_ionized)
@Project.pre(part_4b_job_equilb_NVT_completed_properly)
@Project.pre(part_2c_equilb_NPT_control_file_written)
@Project.post(part_3c_output_equilb_NPT_started)
@Project.post(part_4c_job_equilb_NPT_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get("ncpu"),
        "ngpu": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ngpu", 0
        ),
        "memory": memory_needed,
        "walltime": walltime_gomc_hr,
    }
)
@flow.with_job
@flow.cmd
def run_equilb_NPT_gomc_command(job):
    """Run the gomc equilb_NPT simulation."""
    print("#**********************")
    print("# Started the run_NPT_gomc_command.")
    print("#**********************")

    control_file_name_str = equilb_NPT_control_file_name_str

    print(f"Running simulation job id {job}")
    run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
        str(job.doc.namd_binary_path),
        str(job.doc.equilb_NPT_gomc_binary_file),
        str(ff_info_dict.get(job.sp.forcefield_name).get("ncpu")),
        str(control_file_name_str),
        str(control_file_name_str),
    )

    return run_command


# ******************************************************
# ******************************************************
# equilb NPT or GEMC-NVT - starting the GOMC simulation (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# production run - starting the GOMC simulation (start)
# ******************************************************
# ******************************************************
#@equilibrateSolvent
#This should run the protein system
"""
@Project.pre(part_4c_job_equilb_NPT_completed_properly)
@Project.pre(part_2d_production_control_file_written)
@Project.post(part_3d_output_production_run_started)
@Project.post(part_4d_job_production_run_completed_properly)
@Project.operation.with_directives(
    {
        "np": lambda job: ff_info_dict.get(job.sp.forcefield_name).get("ncpu"),
        "ngpu": lambda job: ff_info_dict.get(job.sp.forcefield_name).get(
            "ngpu", 0
        ),
        "memory": memory_needed,
        "walltime": walltime_gomc_hr,
    }
)
@flow.with_job
def run_production_NPT_gomc_command(job):
    #Run the gomc production_run simulation.
    print("#**********************")
    print("# Started the run_production_NPT_gomc_command function.")
    print("#**********************")

    while (job.doc.cycle < job.doc.num_cycles):
        control_file_name_str = production_control_file_name_str +"_"+str(job.doc.cycle)
    
        print(f"Running simulation job id {job}")
        run_command = "{}/{} +p{} {}.conf > out_{}.dat".format(
            str(job.doc.namd_binary_path),
            str(job.doc.production_ensemble_gomc_binary_file),
            str(ff_info_dict.get(job.sp.forcefield_name).get("ncpu")),
            str(control_file_name_str),
            str(control_file_name_str),
        )

        exec_run_command = subprocess.Popen(
            run_command, shell=True, stderr=subprocess.STDOUT
        )
        os.waitpid(exec_run_command.pid, 0)  # os.WSTOPPED) # 0)

        if part_4d_job_production_run_completed_properly(job):
            job.doc.cycle += 1

"""

def _setup_conf(job, fname, template, data, overwrite=False):

    """Create conf files based on a template and provided data.

    Parameters
    ----------
    fname: str
        Name of the file to be saved out
    template: str, or jinja2.Template
        Either a jinja2.Template or path to a jinja template
    data: dict
        Dictionary storing data matched with the fields available in the template
    overwrite: bool, optional, default=False
        Options to overwrite (or not) existing mdp file of the

    Returns
    -------
    File saved with names defined by fname
    """

    from jinja2 import Template

    if isinstance(template, str):
        with open(template, "r") as f:
            template = Template(f.read())

    if not overwrite:
        if os.path.isfile(fname):
            raise FileExistsError(
                f"{fname} already exists. Set overwrite=True to write out."
            )

    rendered = template.render(data)
    with open(fname, "w") as f:
        f.write(rendered)

    return None


# ******************************************************
# ******************************************************
# production run - starting the GOMC simulation (end)
# ******************************************************
# ******************************************************

# ******************************************************
# ******************************************************
# signac end code (start)
# ******************************************************
# ******************************************************
if __name__ == "__main__":
    pr = Project()
    pr.main()
# ******************************************************
# ******************************************************
# signac end code (end)
# ******************************************************
# ******************************************************
