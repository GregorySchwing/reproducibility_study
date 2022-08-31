"""Setup for signac, signac-flow, signac-dashboard for this study."""
import os
import pathlib
import sys

import flow
import panedr
import unyt as u
from flow.environment import DefaultSlurmEnvironment


from reproducibility_project.src.analysis.equilibration import is_equilibrated
from reproducibility_project.src.utils.forcefields import load_ff
from reproducibility_project.src.proteinffs.forcefields import get_ff_path
from reproducibility_project.src.proteinffs.forcefields import get_wm_path
from reproducibility_project.src.proteinffs.forcefields import get_ff_name
from reproducibility_project.src.proteinffs.forcefields import get_wm_name

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis import hbonds


"""
Pg 3 of 1208351-lindorff-larsen.som.pdf
For each protein we monitored
a set of long-range native Cα–Cα contacts between residues that are separated by at least seven
residues in primary sequence and whose Cα atoms are closer than 10 Å for more than 80% of the
time in the folded state. We then defined a reaction coordinate, Q(t), that monitors the fraction
of these contacts that are formed at each frame interval as the sum:
"""
def contacts_within_cutoff(u, group_a, group_b, radius=10.0, native_threshold=0.80):
    timeseries = []
    contacts_over_time = []
    avg_dist = contacts.distance_array(group_a.positions, group_b.positions)
    avg_dist.fill(0)
    contactmap = contacts.contact_matrix(avg_dist, radius).astype(int)
    contactmap.fill(0)
    # determine which distances <= radius
    #cm = contacts.contact_matrix(dist, radius).astype(int)
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        avg_dist += dist
        # determine which distances <= radius
        this_frame_contactmap = contacts.contact_matrix(dist, radius).astype(int)
        contactmap += this_frame_contactmap
        n_contacts = this_frame_contactmap.sum()
        timeseries.append([ts.frame, n_contacts])

    contactmap = (contactmap.astype('float64') / len(u.trajectory) >= native_threshold)
    avg_dist = avg_dist / len(u.trajectory)
    avg_dist *= contactmap
    avg_dist_norm = (avg_dist - np.min(avg_dist.nonzero())) / (np.max(avg_dist) - np.min(avg_dist.nonzero()))
    avg_dist_norm *= contactmap
    avg_dist_norm[avg_dist_norm == 0.0] = np.nan
    # Make a user-defined colormap.
    cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r","b"])
    cnorm = mcol.Normalize(vmin=np.min(avg_dist.nonzero()),vmax=np.max(avg_dist))

    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)

    #plt.subplot(211)
    
    plt.imshow(avg_dist_norm, cmap=cm1,  interpolation='nearest')
    #plt.subplot(212)
    #plt.imshow(avg_dist_norm, cmap='Greys',  interpolation='nearest')
    plt.colorbar(cpick,label="Native Distance (Å)")
    plt.savefig('native_contacts.png')


    from itertools import product
    atomPairMatrix = np.array(list(product(group_a.atoms.ix_array, group_b.atoms.ix_array)))
    #Flatten 2D contact map rowwise.
    contactmap_flat = contactmap.flatten(order='C')
    from itertools import compress
    # Filter product of c-alpha contacts by flattened contact map
    nativeContacts = list(compress(atomPairMatrix, contactmap_flat))

    # set nonnative contact distances as nan for clarity and checking
    avg_dist[avg_dist == 0.0] = np.nan
    avg_dist_flat = avg_dist.flatten(order='C')
    nativeContactsDistances = list(compress(avg_dist_flat, contactmap_flat))

    #print(nativeContacts)
    plt.clf()
    timeseries_np = np.array(timeseries)
    x, y = zip(*timeseries_np)
    plt.plot(x,y)
    plt.savefig('native_contacts_timeseries.png')

    return (nativeContacts, nativeContactsDistances)

def intermolecular_hbonds(u):
    hbondsa = HydrogenBondAnalysis(universe=u, 
    between=["protein", "protein"],
    d_a_cutoff=5.0,
    d_h_a_angle_cutoff=150,
    update_selections=False)
    hbondsa.run(start=None,
    stop=None,
    step=None,
    verbose=True)    
    print(hbondsa)
    # We see there are 27 hydrogen bonds in total
    print(hbondsa.results.hbonds.shape)
    #Each row contains a unique hydrogen bond. The four columns correspond to:
    #the donor atom index
    #the hydrogen atom index
    #the acceptor atom index
    #the total count
    print(hbondsa.count_by_ids())

    return hbondsa



class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""

    def __init__(self):
        super().__init__()


"""
class Rahman(DefaultPBSEnvironment):

    template = "rahman_gmx.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--walltime",
            type=float,
            default=96,
            help="Walltime for this submission",
        )
 """       

class Grid(DefaultSlurmEnvironment):  # Grid(StandardEnvironment):
    """Subclass of DefaultSlurmEnvironment for WSU's Grid cluster."""

    hostname_pattern = r".*\.grid\.wayne\.edu"
    template = "grid.sh"

import os
import errno
import subprocess
import sys



@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: not j.isfile("solv_ions.gro"))
@Project.pre(lambda j: not j.isfile("topol.top"))
@Project.pre(lambda j: not j.isfile("em.mdp"))
@Project.pre(lambda j: not j.isfile("nvt.mdp"))
@Project.pre(lambda j: not j.isfile("npt_prod.mdp"))
@Project.pre(lambda j: not j.isfile("nvt_prod.mdp"))

@Project.post(lambda j: j.isfile("solv_ions.gro"))
@Project.post(lambda j: j.isfile("topol.top"))
@Project.post(lambda j: j.isfile("em.mdp"))
@Project.post(lambda j: j.isfile("nvt.mdp"))
@Project.post(lambda j: j.isfile("npt_prod.mdp"))
@Project.post(lambda j: j.isfile("nvt_prod.mdp"))
@flow.with_job
def init_job(job):


    """Initialize individual job workspace, including mdp and molecular init files."""
    from reproducibility_project.src.engine_input.gromacs import mdp
    #print("Downloading 6G6J") # 2 Dimers of Myc-Max
    #os.system('pdb_fetch 6G6J > 6G6J.pdb')  # 2.25 Å Resolution
    #print("6G6J Done")
    print("Downloading {}".format(job.sp.pdbid))  # 2.25 Å Resolution  # 1.35 Å Resolution
    os.system("pdb_fetch {} > {}.pdb".format(job.sp.pdbid, job.sp.pdbid))  # 2.25 Å Resolution  # 1.35 Å Resolution
    print("{} Done".format(job.sp.pdbid))  # 2.25 Å Resolution  # 1.35 Å Resolution
    #print("Downloading 6G6L") # 4 Dimers of Myc-Max
    #os.system('pdb_fetch 6G6L > 6G6L.pdb')  # 2.20 Å Resolution
    #print("6G6L Done")

    #1-1-build
    print("Building clean_{} (Stripping Non-protein atoms)".format(job.sp.pdbid))
    os.system("pdb_delhetatm {}.pdb > clean_{}.pdb".format(job.sp.pdbid, job.sp.pdbid))  # 2.25 Å Resolution  # 1.35 Å Resolution)  # 2.25 Å Resolution
    print("clean_{} (Stripping Non-protein atoms) Done".format(job.sp.pdbid))
    # A,C - Myc ; B,D - Max
    # Dimer 1 - A,B
    # Dimer 2 - C,D
    if (job.sp.sim_type == "dimer"):
        job.doc.prot_pdb = "myc_max.pdb"
        job.doc.myc_group = "myc.ndx"
        job.doc.max_group = "max.ndx"

        os.system("pdb_delchain -C,D clean_{}.pdb > {}_single_copy.pdb".format(job.sp.pdbid, job.sp.pdbid))  # 2.25 Å Resolution
        print("myc_max_1.35 (Single Dimer) Done")
    elif (job.sp.sim_type == "myc_mono"):
        job.doc.prot_pdb = "myc_1.35_a_out.pdb"
        os.system("pdb_delchain -B,C,D clean_{}.pdb > {}_single_copy.pdb".format(job.sp.pdbid, job.sp.pdbid))  # 2.25 Å Resolution
        print("myc_1.35 (Single Monomer) Done")
    elif (job.sp.sim_type == "max_mono"):
        job.doc.prot_pdb = "max_1.35_a_out.pdb"
        os.system("pdb_delchain -A,C,D clean_{}.pdb > {}_single_copy.pdb".format(job.sp.pdbid, job.sp.pdbid))  # 2.25 Å Resolution
        print("max_1.35 (Single Monomer) Done")
    else:
        print("Error: bad sim_type")

    job.doc.pH = job.sp.pH
    
    reresCommand = "pdb_reres -1 {}_single_copy.pdb > {}".format(job.sp.pdbid, job.doc.prot_pdb)
    os.system(reresCommand)  # 2.25 Å Resolution

    propPKACommand = "propka3 -f {} --pH {}".format(job.doc.prot_pdb, job.doc.pH)

    linkProtFFCommand = "ln -s {} {}".format(get_ff_path(job.sp.forcefield_name), get_ff_name(job.sp.forcefield_name))
    linkWaterModelCommand = "ln -s {} .".format(get_wm_path(job.sp.forcefield_name))
    pdb2gmxCommand = "gmx pdb2gmx -f {} -o processed.gro <<EOF\n1\n1\nEOF".format(job.doc.prot_pdb)
    pdb2gmxCommandForPLUMED = "gmx pdb2gmx -f {} -o dimer.pdb <<EOF\n1\n1\nEOF".format(job.doc.prot_pdb)

    fixTerminiForPLUMEDCommand = "find dimer.pdb -type f -exec sed -i 's/OC1/O  /g' {} \;"
    #"find dimer_bad_termini.pdb -type f -exec sed -i 's/OC2/OXT/g' {} \;"

    print(linkProtFFCommand)
    print(linkWaterModelCommand)
    print(reresCommand)
    print(pdb2gmxCommand)
    print(pdb2gmxCommandForPLUMED)
    print(fixTerminiForPLUMEDCommand)

    os.system(linkProtFFCommand)  # 2.25 Å Resolution
    os.system(linkWaterModelCommand)  # 2.25 Å Resolution
    ####Make a topology file using structure and force field for simulation. Make sure to have a structure file of a protein (e.g., histatin5.pdb) and a force field directory if one is using a different force field other than the available in the compiled version of the gromacs. pdb2gmx asks to choose a force field and water model. In this example, it will choose the force field and water model listed in option 1. Check and make sure.
    os.system(pdb2gmxCommand)  # 2.25 Å Resolution
    os.system(pdb2gmxCommandForPLUMED)  # 2.25 Å Resolution
    os.system(fixTerminiForPLUMEDCommand)  # 2.25 Å Resolution


    ####Prepare a simulaton box. For IDP, box dimension need to be large enough to prevent any periodic image interaction.
    editconfCommand = "gmx editconf -f processed.gro -o newbox.gro -c -d 2.0 -bt cubic"
    os.system(editconfCommand)  # 2.25 Å Resolution


    ####Make sure to have water structure file (e.g., tip4p2005.gro) in the working directory.
    ####Solvating a simulation box.
    solvateCommand = "gmx solvate -cp newbox.gro -cs {} -o solv.gro -p topol.top".format(get_wm_name(job.sp.forcefield_name))
    os.system(solvateCommand)  # 2.25 Å Resolution

    ####Adding counter ions to neutralize the box. Replace "SOL" while adding ions.
    grommpPreGenionCommand = "gmx grompp -f {}/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1".format(os.path.dirname(os.path.abspath(mdp.__file__)))
    if (job.sp.salt_conc == "None"):
        genionCommand = "echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname {} -pq {} -nname {} -nq {} -neutral -seed {}".format(job.sp.cat_name, job.sp.cat_val, job.sp.an_name, job.sp.an_val, job.sp.genion_seed)
    else:
        genionCommand = "echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname {} -pq {} -nname {} -nq {} -neutral -conc {} -seed {}".format(job.sp.cat_name, job.sp.cat_val, job.sp.an_name, job.sp.an_val, job.sp.salt_conc, job.sp.genion_seed)

    os.system(grommpPreGenionCommand)  # 2.25 Å Resolution
    os.system(genionCommand)  # 2.25 Å Resolution

    pressure = job.sp.pressure * u.kPa
    mdp_abs_path = os.path.dirname(os.path.abspath(mdp.__file__))

    
    if (job.sp.sim_type == "dimer"):
        makeNDX = "gmx make_ndx -f solv_ions.gro -o index.ndx <<EOF\nr 1-89\nname 19 myc\nr 90-165\nname 20 max\nquit\nEOF"
        #makeNDX = "gmx make_ndx -f {} -o index.ndx <<EOF\nchain A\nchain B\nquit\nEOF".format(job.doc.prot_pdb)

        os.system(makeNDX)  # 2.25 Å Resolution

    elif (job.sp.sim_type == "myc_mono"):
        makeNDX = "gmx make_ndx -f {} -o index.ndx <<EOF\nchain A\nquit\nEOF".format(job.doc.prot_pdb)
        os.system(makeNDX)  # 2.25 Å Resolution

    elif (job.sp.sim_type == "max_mono"):
        makeNDX = "gmx make_ndx -f {} -o index.ndx <<EOF\nchain B\nquit\nEOF".format(job.doc.prot_pdb)
        os.system(makeNDX)  # 2.25 Å Resolution
    

    mdps = {
        "em": {
            "fname": "em.mdp",
            "template": f"{mdp_abs_path}/em_template.mdp.jinja",
            "water-template": f"{mdp_abs_path}/em_template_water.mdp.jinja",
            "data": {
                "temp": job.sp.temperature,
                "replica": job.sp.replica,
            },
        },
        "nvt_eq": {
            "fname": "nvt_eq.mdp",
            "template": f"{mdp_abs_path}/nvt_template.mdp.jinja",
            "water-template": f"{mdp_abs_path}/nvt_template_water.mdp.jinja",
            "data": {
                #"nsteps": 2500000,
                "nsteps": 250,
                "dt": 0.002,
                "temp": job.sp.temperature,
            },
        },
        "npt_eq": {
            "fname": "npt_eq.mdp",
            "template": f"{mdp_abs_path}/npt_template.mdp.jinja",
            "water-template": f"{mdp_abs_path}/npt_template_water.mdp.jinja",
            "data": {
                #"nsteps": 5000000,
                "nsteps": 500,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "refp": pressure.to_value("bar"),
            },
        },
        "npt_prod": {
            "fname": "npt_prod.mdp",
            "template": f"{mdp_abs_path}/npt_prod_template.mdp.jinja",
            "water-template": f"{mdp_abs_path}/npt_template_water.mdp.jinja",
            "data": {
                #"nsteps": 5000000,
                "nsteps": 500,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "refp": pressure.to_value("bar"),
            },
        },
        "nvt_prod": {
            "fname": "nvt_prod.mdp",
            "template": f"{mdp_abs_path}/nvt_template.mdp.jinja",
            "water-template": f"{mdp_abs_path}/nvt_template_water.mdp.jinja",
            "data": {
                #"nsteps": 5000000,
                "nsteps": 500,
                "dt": 0.001,
                "temp": job.sp.temperature,
            },
        },
        "npt_plumed_metadynamics": {
            "fname": "npt_plumed_metadynamics.mdp",
            "template": f"{mdp_abs_path}/npt_plumed_metadynamics_template.mdp.jinja",
            "water-template": f"{mdp_abs_path}/npt_plumed_metadynamics_template_water.mdp.jinja",
            "data": {
                #"nsteps": 5000000,
                "nsteps": 500,
                "dt": 0.001,
                "temp": job.sp.temperature,
                "refp": pressure.to_value("bar"),
            },
        },
    }

    for op, mdp in mdps.items():
        _setup_mdp(
            fname=mdp["fname"],
            template=mdp["template"],
            data=mdp["data"],
            overwrite=True,
        )


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("solv_ions.gro"))
@Project.pre(lambda j: j.isfile("topol.top"))
@Project.post(lambda j: j.isfile("em.gro"))
@flow.with_job
@flow.cmd
def gmx_em(job):
    """Run GROMACS grompp for the energy minimization step."""
    em_mdp_path = "em.mdp"
    grompp = f"gmx grompp -f {em_mdp_path} -o em.tpr -c solv_ions.gro -p topol.top --maxwarn 1"
    mdrun = _mdrun_str("em")
    return f"{grompp} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("em.gro"))
@Project.post(lambda j: j.isfile("nvt_eq.gro"))
@flow.with_job
@flow.cmd
def gmx_nvt_eq(job):
    """Run GROMACS grompp for the nvt_eq step."""
    nvt_eq_mdp_path = "nvt_eq.mdp"
    grompp = f"gmx grompp -f {nvt_eq_mdp_path} -o nvt_eq.tpr -c em.gro -r em.gro -p topol.top --maxwarn 1"
    mdrun = _mdrun_str("nvt_eq")
    return f"{grompp} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt_eq.gro"))
@Project.post(lambda j: j.isfile("npt_eq.gro"))
@flow.with_job
@flow.cmd
def gmx_npt_eq(job):
    """Run GROMACS grompp for the npt step."""
    npt_eq_mdp_path = "npt_eq.mdp"
    grompp = f"gmx grompp -f {npt_eq_mdp_path} -o npt_eq.tpr -c nvt_eq.gro -r nvt_eq.gro -p topol.top --maxwarn 1"
    mdrun = _mdrun_str("npt_eq")
    return f"{grompp} && {mdrun}"

@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_eq.gro"))
@Project.post(lambda j: j.isfile("npt_prod.gro"))
@flow.with_job
@flow.cmd
def gmx_npt_prod(job):
    """Run GROMACS grompp for the npt step."""
    npt_prod_mdp_path = "npt_prod.mdp"
    grompp = f"gmx grompp -f {npt_prod_mdp_path} -o npt_prod.tpr -c npt_eq.gro -p topol.top --maxwarn 1"
    mdrun = _mdrun_str("npt_prod")
    return f"{grompp} && {mdrun}"

"""
This method currently requires a trr of the structure, since native and crystal
likely differ.  Though https://pubs.acs.org/action/showCitFormats?doi=10.1021/jacs.0c03217&ref=pdf
used the crystal to define the natives.
while the more robust https://www.jstor.org/stable/41351322
defined natives as < 10 A on average of 80% of the conformations. 
"""

@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.trr"))
@Project.pre(lambda j: j.isfile("npt_prod.tpr"))
@Project.pre(lambda j: "native_contact_list" not in j.doc)
#@Project.pre(lambda j: equil_status(j, "npt_prod", "Volume"))
@Project.post(lambda j: "native_contact_list" in j.doc)
@flow.with_job
def determine_native_dimer_contacts(job):
    u = mda.Universe(job.fn("npt_prod.tpr"), job.fn("npt_prod.trr"))

    sel_myc = "resid 1:89 and name CA"
    sel_max = "resid 90:165 and name CA"
    # reference groups (first frame of the trajectory, but you could also use a
    # separate PDB, eg crystal structure)
    myc = u.select_atoms(sel_myc)
    max = u.select_atoms(sel_max)

    cts = contacts_within_cutoff(u, myc, max)
    job.doc.native_contact_list = cts[0]
    job.doc.native_contact_average_distances = cts[1]


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.pre(lambda j: j.isfile("plumed.dat"))
@Project.post(lambda j: j.isfile("npt_plumed_metadynamics.gro"))

@flow.with_job
@flow.cmd
def gmx_npt_plumed_metadynamics(job):
    """Run GROMACS grompp for the npt step."""
    npt_plumed_metadynamics_mdp_path = "npt_plumed_metadynamics.mdp"
    grompp = f"gmx grompp -f {npt_plumed_metadynamics_mdp_path} -o npt_plumed_metadynamics.tpr -c npt_prod.gro -p topol.top --maxwarn 1"
    mdrun = _mdrun_plumed_str("npt_plumed_metadynamics")
    return f"{grompp} && {mdrun}"

"""

@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.pre(
    lambda j: not equil_status(j, "npt_prod", "Potential")
    or not equil_status(j, "npt_prod", "Volume")
)
@Project.post(lambda j: equil_status(j, "npt_prod", "Potential"))
@Project.post(lambda j: equil_status(j, "npt_prod", "Volume"))
@flow.with_job
@flow.cmd
def extend_gmx_npt_prod(job):
    #Run GROMACS grompp for the npt step.
    # Extend the npt run by 1000 ps (1 ns)
    extend = "gmx convert-tpr -s npt_prod.tpr -extend 1000 -o npt_prod.tpr"
    mdrun = _mdrun_str("npt_prod")
    return f"{extend} && {mdrun}"

@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.pre(
    lambda j: not equil_status(j, "npt_prod", "Potential")
    or not equil_status(j, "npt_prod", "Volume")
)
@Project.post(lambda j: equil_status(j, "npt_prod", "Potential"))
@Project.post(lambda j: equil_status(j, "npt_prod", "Volume"))
@flow.with_job
@flow.cmd
def extend_gmx_npt_prod(job):
    #Run GROMACS grompp for the npt step.
    # Extend the npt run by 1000 ps (1 ns)
    extend = "gmx convert-tpr -s npt_prod.tpr -extend 1000 -o npt_prod.tpr"
    mdrun = _mdrun_str("npt_prod")
    return f"{extend} && {mdrun}"



@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.post(lambda j: j.isfile("nvt_prod.gro"))
@flow.with_job
@flow.cmd
def gmx_nvt_prod(job):
    #Run GROMACS grompp for the nvt step.
    npt_mdp_path = "npt_prod.mdp"
    grompp = f"gmx grompp -f {npt_mdp_path} -o nvt_prod.tpr -c npt_prod.gro -p topol.top --maxwarn 1"
    mdrun = _mdrun_str("nvt_prod")
    return f"{grompp} && {mdrun}"


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt_prod.gro"))
@Project.pre(
    lambda j: not equil_status(j, "nvt_prod", "Potential")
    or not equil_status(j, "nvt_prod", "Pressure")
)
@Project.post(lambda j: equil_status(j, "nvt_prod", "Potential"))
@Project.post(lambda j: equil_status(j, "nvt_prod", "Pressure"))
@flow.with_job
@flow.cmd
def extend_gmx_nvt_prod(job):
    # Run GROMACS grompp for the nvt step.
    # Extend the nvt run by 1000 ps (1 ns)
    extend = "gmx convert-tpr -s nvt_prod.tpr -extend 1000 -o nvt_prod.tpr"
    mdrun = _mdrun_str("nvt_prod")
    return f"{extend} && {mdrun}"
"""



@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.pre(lambda j: j.isfile("npt_prod.tpr"))
@Project.pre(lambda j: "native_contact_list" in j.doc)
@Project.pre(lambda j: "native_contact_average_distances" in j.doc)
@Project.pre(lambda j: not j.isfile("plumed.dat"))
@Project.post(lambda j: j.isfile("plumed.dat"))

@flow.with_job
def create_plumed_file(job):
    from reproducibility_project.src.engine_input.gromacs import mdp
    from string import Template
    t = Template("\tATOMS$id=$atom1,$atom2 SWITCH$id={$_SWITCHFUNCTION R_0=0.01 BETA=$beta LAMBDA=$_lambda REF=$dist}  WEIGHT$id=$weight\n")
    contact_map_string = ""
    counter = 1
    weight = 1/len(job.doc.native_contact_average_distances)
    _SWITCHFUNCTION="Q"
    # Refer to Eq (1) of https://pubs.acs.org/action/showCitFormats?doi=10.1021/jacs.0c03217&ref=pdf
    beta=10.0
    _lambda=1.0
    for f, b in zip(job.doc.native_contact_list, job.doc.native_contact_average_distances):
        string_template = t.substitute(id=counter,atom1=f[0],atom2=f[1], _SWITCHFUNCTION=_SWITCHFUNCTION, beta=beta, _lambda=_lambda, dist=b, weight=weight)
        contact_map_string += string_template
        counter = counter + 1
    #print(contact_map_string)
    mdp_abs_path = os.path.dirname(os.path.abspath(mdp.__file__))
    mdps = {
        "plumed": {
            "fname": "plumed.dat",
            "template": f"{mdp_abs_path}/test.dat.jinja",
            "water-template": f"{mdp_abs_path}/plumed.dat.mdp.jinja",
            "data": {
                "myc_res": "1-89",
                "max_res": "90-165",
                "cmstring": contact_map_string,
                #"temp": job.sp.temperature,
                #"refp": pressure.to_value("bar"),
            },
        },
    }
    for op, mdp in mdps.items():
        _setup_mdp(
            fname=mdp["fname"],
            template=mdp["template"],
            data=mdp["data"],
            overwrite=True,
        )
    

@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_plumed_metadynamics.gro"))
#@Project.pre(lambda j: equil_status(j, "npt_prod", "Potential"))
#@Project.pre(lambda j: equil_status(j, "npt_prod", "Volume"))
@Project.post(lambda j: j.isfile("log-npt.txt"))
@Project.post(lambda j: j.isfile("trajectory-npt.gsd"))
@flow.with_job
def sample_protein_properties(job):
    """Sample properties of interest from npt edr."""
    import mdtraj
    import pandas as pd

    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
        sample_job,
    )

    p = pathlib.Path(job.workspace())
    data = panedr.edr_to_df(f"{str(p.absolute())}/npt.edr")
    # Properties of interest
    poi = {
        "Potential": "potential_energy",
        "Kinetic En.": "kinetic_energy",
        "Pressure": "pressure",
        "Temperature": "temperature",
        "Density": "density",
        "Volume": "volume",
    }

    tmp_df = pd.DataFrame()
    for idx, row in data.iterrows():
        tmp_df = tmp_df.append(row[list(poi.keys())])
    tmp_df.rename(poi, axis=1, inplace=True)
    tmp_df.insert(
        column="time_steps",
        value=[10000 * i for i in range(len(tmp_df))],
        loc=1,
    )
    tmp_df.to_csv("log-npt.txt", index=False, sep=" ")
    for prop in poi:
        sample_job(job, filename="log-npt.txt", variable=poi[prop])
        get_subsampled_values(
            job,
            property=poi[prop],
            property_filename="log-npt.txt",
            ensemble="npt",
        )

    # Convert trr file to gsd with mdtraj
    traj = mdtraj.load("npt_prod.trr", top="npt_prod.gro")
    traj.save("trajectory-npt.gsd")


"""
@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("npt_prod.gro"))
@Project.pre(lambda j: equil_status(j, "npt_prod", "Potential"))
@Project.pre(lambda j: equil_status(j, "npt_prod", "Volume"))
@Project.post(lambda j: j.isfile("log-npt.txt"))
@Project.post(lambda j: j.isfile("trajectory-npt.gsd"))
@flow.with_job
def sample_npt_properties(job):
    #Sample properties of interest from npt edr.
    import mdtraj
    import pandas as pd

    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
        sample_job,
    )

    p = pathlib.Path(job.workspace())
    data = panedr.edr_to_df(f"{str(p.absolute())}/npt.edr")
    # Properties of interest
    poi = {
        "Potential": "potential_energy",
        "Kinetic En.": "kinetic_energy",
        "Pressure": "pressure",
        "Temperature": "temperature",
        "Density": "density",
        "Volume": "volume",
    }

    tmp_df = pd.DataFrame()
    for idx, row in data.iterrows():
        tmp_df = tmp_df.append(row[list(poi.keys())])
    tmp_df.rename(poi, axis=1, inplace=True)
    tmp_df.insert(
        column="time_steps",
        value=[10000 * i for i in range(len(tmp_df))],
        loc=1,
    )
    tmp_df.to_csv("log-npt.txt", index=False, sep=" ")
    for prop in poi:
        sample_job(job, filename="log-npt.txt", variable=poi[prop])
        get_subsampled_values(
            job,
            property=poi[prop],
            property_filename="log-npt.txt",
            ensemble="npt",
        )

    # Convert trr file to gsd with mdtraj
    traj = mdtraj.load("npt_prod.trr", top="npt_prod.gro")
    traj.save("trajectory-npt.gsd")


@Project.operation
@Project.pre(lambda j: j.sp.engine == "gromacs")
@Project.pre(lambda j: j.isfile("nvt_prod.gro"))
@Project.pre(lambda j: equil_status(j, "nvt_prod", "Potential"))
@Project.pre(lambda j: equil_status(j, "nvt_prod", "Pressure"))
@Project.post(lambda j: j.isfile("log-nvt.txt"))
@Project.post(lambda j: j.isfile("trajectory-nvt.gsd"))
@flow.with_job
def sample_nvt_properties(job):
    #Sample properties of interest from nvt edr.
    import mdtraj
    import pandas as pd

    from reproducibility_project.src.analysis.sampler import (
        get_subsampled_values,
        sample_job,
    )

    p = pathlib.Path(job.workspace())
    data = panedr.edr_to_df(f"{str(p.absolute())}/nvt_prod.edr")
    # Properties of interest
    poi = {
        "Time": "time",
        "Potential": "potential_energy",
        "Kinetic En.": "kinetic_energy",
        "Pressure": "pressure",
        "Temperature": "temperature",
    }

    tmp_df = pd.DataFrame()
    for idx, row in data.iterrows():
        tmp_df = tmp_df.append(row[list(poi.keys())])
    tmp_df.rename(poi, axis=1, inplace=True)
    tmp_df.insert(
        column="time_steps",
        value=[10000 * i for i in range(len(tmp_df))],
        loc=1,
    )
    tmp_df.to_csv("log-nvt.txt", index=False, sep=" ")
    for prop in poi:
        sample_job(job, filename="log-nvt.txt", variable=poi[prop])
        get_subsampled_values(
            job,
            property=poi[prop],
            property_filename="log-nvt.txt",
            ensemble="nvt",
        )

    # Convert trr file to gsd with mdtraj
    traj = mdtraj.load("nvt_prod.trr", top="nvt_prod.gro")
    traj.save("trajectory-nvt.gsd")
"""

def _mdrun_str(op):
    """Output an mdrun string for arbitrary operation."""
    #PLUMED 2021.5 is Slow
    #msg = f"/usr/local/gromacs_2022/bin/gmx mdrun -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt -nt 16"
    msg = f"/usr/local/gromacs_2022/bin/gmx mdrun -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt -nt 16"
    return msg

def _mdrun_plumed_str(op):
    """Output an mdrun string for arbitrary operation."""
    msg = f"gmx mdrun -plumed plumed.dat -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt -nt 16"
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


def equil_status(job, op, att):
    """Check equilibration status of specific attributes of specific operation."""
    p = pathlib.Path(job.workspace())
    if not job.isfile(f"{op}.edr"):
        return False
    else:
        data = panedr.edr_to_df(f"{str(p.absolute())}/{op}.edr")
        return is_equilibrated(data[att])[0]


if __name__ == "__main__":
    pr = Project()
    pr.main()
