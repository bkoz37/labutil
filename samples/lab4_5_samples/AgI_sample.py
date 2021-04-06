import numpy, os
import matplotlib.pyplot as plt
from labutil.plugins.lammps import lammps_run, parse_lammps_rdf, parse_lammps_thermo, get_rdf
from labutil.objects import File, Struc, Dir, ase2struc, ClassicalPotential
from ase.spacegroup import crystal
from ase.build import make_supercell
from ase import Atoms

intemplate = """
    # ---------- Initialize simulation ---------------------
    units metal
    atom_style atomic
    dimension  3
    boundary   p p p
    newton off
    read_data $DATAINPUT 

    pair_style mff
    pair_coeff * * /home/bond/Work/Lab4/AgI_FF.txt 47 53 yes yes

    velocity all create $TEMPERATURE 126342 dist gaussian rot yes mom yes

    group silvers type 1
    group iodines type 2

    # ---------- Describe computed properties------------------
    compute msdAg silvers msd com yes
    compute msdI iodines msd com yes
    compute rdfAg all rdf 1000 1 1
    compute rdfI all rdf 1000 2 2
    compute rdfAgI all rdf 1000 1 2

    variable rdfAgFile string "$RDFFILE.Ag"
    variable rdfIFile string "$RDFFILE.I"
    variable rdfAgIFile string "$RDFFILE.AgI"

    thermo_style custom step temp etotal press density c_msdAg[4] c_msdI[4]
    thermo $TOUTPUT 

    # record rdf
    fix 1 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfAg[*] file ${rdfAgFile} mode vector
    fix 2 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfI[*] file ${rdfIFile} mode vector
    fix 3 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfAgI[*] file ${rdfAgIFile} mode vector

    # ---------- Specify ensemble  ---------------------
    fix 4 all npt temp $TEMPERATURE $TEMPERATURE $TDAMP tri 0.0 0.0 1.0

    # --------- Run -------------
    timestep $TIMESTEP
    run $NSTEPS
    unfix 4
"""

def make_struc(size):
    """
    Creates the crystal structure using ASE.
    :param size: supercell multiplier
    :return: structure object converted from ase
    """
    alat = 5.1
    lattice = alat * numpy.identity(3)
    symbols = ['Ag', 'I', 'Ag', 'I']
    sc_positions = [[1/2, 0, 1/4], [0, 0, 0], [1, 1/2, 3/4], [1/2, 1/2, 1/2]]   
    unitcell = Atoms(symbols=symbols, scaled_positions=sc_positions, cell=lattice)
    multiplier = numpy.identity(3) * size
    supercell = make_supercell(unitcell, multiplier)
    structure = Struc(ase2struc(supercell))
    return structure


def compute_AgI_dynamics(size, timestep, nsteps, temperature, ncpu):
    """
    Make an input template and select potential and structure, and input parameters.
    Return a pair of output file and RDF file written to the runpath directory.
    """

    potential = ClassicalPotential(ptype='eam', element='Al', name='Al_zhou.eam.alloy')
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab4/Problem2", "temp_" + str(temperature)))
    struc = make_struc(size=size)
    inparam = {
        'TEMPERATURE': temperature,
        'NSTEPS': nsteps,
        'TIMESTEP': timestep,
        'TOUTPUT': 100,                 # how often to write thermo output
        'TDAMP': 50 * timestep,       # thermostat damping time scale
        'RDFFRAME': int(nsteps / 4),   # frames for radial distribution function
    }
    outfile = lammps_run(struc=struc, runpath=runpath, potential=potential,
                                  intemplate=intemplate, inparam=inparam, ncpu=ncpu, triclinic=True)
    output = parse_lammps_thermo(outfile=outfile)

    rdfAgFile = File(path=os.path.join(runpath.path, 'lammps.rdf.Ag')) 
    rdfIFile = File(path=os.path.join(runpath.path, 'lammps.rdf.I')) 
    rdfAgIFile = File(path=os.path.join(runpath.path, 'lammps.rdf.AgI')) 

    rdfsAg = parse_lammps_rdf(rdffile=rdfAgFile)
    rdfsI = parse_lammps_rdf(rdffile=rdfIFile)
    rdfsAgI = parse_lammps_rdf(rdffile=rdfAgIFile)

    return output, rdfsAg, rdfsI, rdfsAgI


def md_run():
    output, rdfsAg, rdfsI, rdfsAgI = compute_AgI_dynamics(size=1, timestep=0.001, nsteps=1000, temperature=300, ncpu=1)
    [simtime, temp, etotal, press, dens, msdAg, msdI] = output
    ## ------- plot output properties
    plt.plot(simtime, temp)
    plt.show()

    # ----- plot radial distribution functions
    for rdf in rdfsAgI:
        plt.plot(rdf[0], rdf[1])
    plt.show()


if __name__ == '__main__':
    # put here the function that you actually want to run
    md_run()
