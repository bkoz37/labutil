import numpy, os
import matplotlib.pyplot as plt
from labutil.plugins.lammps import lammps_run, parse_lammps_rdf, parse_lammps_thermo, get_rdf
from labutil.objects import File, Struc, Dir, ase2struc, ClassicalPotential
from ase.spacegroup import crystal
from ase.build import make_supercell

intemplate = """
    # ---------- Initialize simulation ---------------------
    units metal
    atom_style atomic
    dimension  3
    boundary   p p p
    newton off
    read_data /home/bond/Work/Lab4/AgI.data 

    pair_style mgpf
    pair_coeff * * /home/bond/Work/Lab4/AgI_FF.txt 47 53 yes yes

    velocity all create $TEMPERATURE 126342 dist gaussian rot yes mom yes

    group silvers type 1
    group iodines type 2

    # ---------- Describe computed properties------------------
    compute msdAg silvers msd com yes
    compute msdI iodines msd com yes
    compute rdfAg all rdf 1000 1 1
    compute rdfI all rdf 1000 2 2

    variable rdfAgFile string "$RDFFILE.Ag"
    variable rdfIFile string "RDFFILE.I"

    thermo_style custom step temp etotal press density c_msdAg[4] c_msdI[4]
    thermo $TOUTPUT 

    # record rdf
    fix 1 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfAg[*] file ${rdfAgFile} mode vector
    fix 1 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfI[*] file ${rdfIFile} mode vector

    # ---------- Specify ensemble  ---------------------
    fix 4 all npt temp $TEMPERATURE $TEMPERATURE $TDAMP tri 0.0 0.0 1.0

    # --------- Run -------------
    run $NSTEPS
    unfix 4
"""

def make_struc(size):
    """
    Creates the crystal structure using ASE.
    :param size: supercell multiplier
    :return: structure object converted from ase
    """
    alat = 4.10
    unitcell = crystal('Al', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
    multiplier = numpy.identity(3) * size
    supercell = make_supercell(unitcell, multiplier)
    structure = Struc(ase2struc(supercell))
    return structure

def compute_AgI_dynamics(timestep, nsteps, temperature, ncpu):
    """
    Make an input template and select potential and structure, and input parameters.
    Return a pair of output file and RDF file written to the runpath directory.
    """

    potential = ClassicalPotential(ptype='eam', element='Al', name='Al_zhou.eam.alloy')
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab4/Problem2", "temp_" + str(temperature)))
    struc = make_struc(size=1)
    inparam = {
        'TEMPERATURE': temperature,
        'NSTEPS': nsteps,
        'TIMESTEP': timestep,
        'TOUTPUT': 100,                 # how often to write thermo output
        'TDAMP': 50 * timestep,       # thermostat damping time scale
        'RDFFRAME': int(nsteps / 4),   # frames for radial distribution function
    }
    outfile = lammps_run(struc=struc, runpath=runpath, potential=potential,
                                  intemplate=intemplate, inparam=inparam, ncpu=ncpu)
    output = parse_lammps_thermo(outfile=outfile)
    rdfAgFile = File(path=os.path.join(runpath.path, 'lammps.rdf.Ag')) 
    rdfIFile = File(path=os.path.join(runpath.path, 'lammps.rdf.I')) 
    rdfsAg = parse_lammps_rdf(rdffile=rdfAgFile)
    rdfsI = parse_lammps_rdf(rdffile=rdfIFile)

    return output, rdfsAg, rdfsI


def md_run():
    output, rdfsAg, rdfsI = compute_AgI_dynamics(timestep=0.001, nsteps=1000, temperature=300, ncpu=1)
    [simtime, temp, etotal, press, dens, msdAg, msdI] = output
    ## ------- plot output properties
    #plt.plot(simtime, temp)
    #plt.show()
    plt.plot(simtime, press)
    plt.show()

    # ----- plot radial distribution functions
    for rdf in rdfsAg:
        plt.plot(rdf[0], rdf[1])
    plt.show()


if __name__ == '__main__':
    # put here the function that you actually want to run
    md_run()
