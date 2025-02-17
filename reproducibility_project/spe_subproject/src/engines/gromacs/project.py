"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow
import pandas as pd
import panedr
import unyt as u
from flow.environment import DefaultPBSEnvironment

from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.utils.forcefields import load_ff


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.post(lambda j: j.isfile("init.gro"))
@Project.post(lambda j: j.isfile("init.top"))
@Project.post(lambda j: j.isfile("nvt.mdp"))
@Project.post(lambda j: j.isfile("npt.mdp"))
@flow.with_job
def LoadSystemSnapShot(job):
    """Create initial configurations of the system statepoint."""
    import mbuild as mb

    from reproducibility_project.spe_subproject.src.engine_input.gromacs import (
        mdp,
    )

    pr = Project()
    snapshot_directory = (
        pathlib.Path(pr.root_directory()) / "src" / "system_snapshots"
    )
    molecule = job.sp.molecule
    molecule_filename = molecule + ".json"
    box = mb.load(str(snapshot_directory / molecule_filename))
    box.save(filename="init.gro", precision=8, overwrite=True)
    # Apply forcefield and write out engine input files
    # __________________________________________________
    ff = load_ff(job.sp.forcefield_name)
    param_system = ff.apply(box)
    param_system.save(
        "init.top",
        overwrite=True,
    )

    # Modify mdp files according to job statepoint parameters
    cutoff_styles = {"hard": "None", "shift": "Potential-shift"}
    lrcs = {"None": "no", "energy_pressure": "EnerPres"}

    pressure = job.sp.pressure * u.kPa
    mdp_abs_path = os.path.dirname(os.path.abspath(mdp.__file__))
    mdps = {
        "npt": {
            "fname": "npt.mdp",
            "template": f"{mdp_abs_path}/npt_template.mdp.jinja",
            "water-template": f"{mdp_abs_path}/npt_template_water.mdp.jinja",
            "data": {
                "nsteps": 1,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "refp": pressure.to_value("bar"),
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "lrc": lrcs[job.sp.long_range_correction],
            },
        },
        "nvt": {
            "fname": "nvt.mdp",
            "template": f"{mdp_abs_path}/nvt_template.mdp.jinja",
            "water-template": f"{mdp_abs_path}/nvt_template_water.mdp.jinja",
            "data": {
                "nsteps": 1,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "r_cut": job.sp.r_cut,
                "cutoff_style": cutoff_styles[job.sp.cutoff_style],
                "lrc": lrcs[job.sp.long_range_correction],
            },
        },
    }

    for op, mdp in mdps.items():
        if job.sp.molecule == "waterSPCE":
            _setup_mdp(
                fname=mdp["fname"],
                template=mdp["water-template"],
                data=mdp["data"],
                overwrite=True,
            )
        else:
            _setup_mdp(
                fname=mdp["fname"],
                template=mdp["template"],
                data=mdp["data"],
                overwrite=True,
            )


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("init.gro"))
@Project.pre(lambda j: j.isfile("init.top"))
@Project.pre(lambda j: j.isfile("nvt.mdp"))
@Project.pre(lambda j: j.isfile("npt.mdp"))
@Project.post(lambda j: j.isfile("nvt.edr"))
@flow.with_job
@flow.cmd
def CalculateEnergy(job):
    """Load onto a cluster and output the point energy for the snapshot."""
    nvt_mdp_path = "nvt.mdp"
    grompp = f"gmx grompp -f {nvt_mdp_path} -o nvt.tpr -c init.gro -p init.top --maxwarn 1"
    mdrun = _mdrun_str("nvt")
    return f"{grompp} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt.edr"))
@Project.post(lambda j: j.isfile("log-spe.txt"))
@flow.with_job
def FormatTextFile(job):
    """Take the output from the simulation engine and convert it to log-spe.txt for data comparisons.

    See README.md for spe_subproject for formatting information.
    """
    import mdtraj
    import panedr

    p = pathlib.Path(job.workspace())
    data = panedr.edr_to_df(f"{str(p.absolute())}/nvt.edr").loc[0.0]

    to_drop = [
        "Vir-XX",
        "Vir-XY",
        "Vir-XZ",
        "Vir-YX",
        "Vir-YY",
        "Vir-YZ",
        "Vir-ZX",
        "Vir-ZY",
        "Vir-ZZ",
        "Pres-XX",
        "Pres-XY",
        "Pres-XZ",
        "Pres-YX",
        "Pres-YY",
        "Pres-YZ",
        "Pres-ZX",
        "Pres-ZY",
        "Pres-ZZ",
        "#Surf*SurfTen",
        "T-System",
        "Conserved En.",
        "Temperature",
        "Pres. DC",
        "Pressure",
        "Total Energy",
    ]
    for key in to_drop:
        data.pop(key)

    spe = {
        "potential_energy": data.get("Potential"),
        "vdw_energy": data.get("LJ (SR)"),
        "tail_energy": data.get("Disper. corr."),
        "coul_energy": data.get("Coulomb (SR)"),
        "kspace_energy": data.get("Coul. recip."),
        "pair_energy": _get_gmx_energy_pair(dict(data)),
        "bonds_energy": data.get("Bond"),
        "angles_energy": data.get("Angle"),
        "dihedrals_energy": _get_gmx_energy_torsion(dict(data)),
    }
    spe_df = pd.DataFrame(spe, index=[0])
    spe_df.to_csv("log-spe.txt", header=True, index=False, sep=",")


"""
The below methods are adapted from
https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/drivers/gromacs.py
"""


def _get_gmx_energy_pair(gmx_energies):
    gmx_pairs = 0.0
    for key in ["LJ-14", "Coulomb-14"]:
        try:
            gmx_pairs += gmx_energies[key]
        except KeyError:
            pass
    return gmx_pairs


def _get_gmx_energy_torsion(gmx_energies):
    """Canonicalize torsion energies from a set of GROMACS energies."""
    gmx_torsion = 0.0
    for key in ["Torsion", "Ryckaert-Bell.", "Proper Dih."]:
        try:
            gmx_torsion += gmx_energies[key]
        except KeyError:
            pass

    return gmx_torsion


def _mdrun_str(op):
    """Output an mdrun string for arbitrary operation."""
    msg = f"gmx mdrun -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt -nt 16"
    return msg


def _setup_mdp(fname, template, data, overwrite=False):
    """Create mdp files based on a template and provided data.

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


if __name__ == "__main__":
    pr = Project()
    pr.main()
