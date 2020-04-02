import numpy, os
import matplotlib.pyplot as plt
from labutil.plugins.lammps import lammps_run, parse_lammps_rdf, parse_lammps_thermo, get_rdf
from labutil.objects import Struc, Dir, ase2struc, ClassicalPotential
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

    pair_style mff
    pair_coeff * * /home/bond/Work/Lab4/AgI_Three.txt 47 53 yes yes

    velocity all create $TEMPERATURE 126342 dist gaussian rot yes mom yes

    thermo_style one
    thermo $TOUTPUT 

    #dump 1 all custom 100000 ./AgI.dump id type x y z fx fy fz vx vy vz
    #dump_modify 1 sort id

    group silvers type 1
    group iodines type 2

    # ---------- Describe computed properties------------------
    compute msdall all msd com yes
    compute msdAg silvers msd com yes
    compute msdI iodines msd com yes
    compute rdfall all rdf 1000
    compute rdfpairs all rdf 1000 1 1 1 2 2 2

    variable v equal vol
    variable t equal temp

    # record rdf and msd
    fix 1 all ave/time 1 $RDFFRAME $RDFFRAME c_rdfall[*] c_rdfpairs[*] file $RDFFILE mode vector
    fix 2 all ave/time 10000 1 10000 c_msdall[*] c_msdAg[*] c_msdI[*] file AgI.msd mode vector

    # record rolling average of volume and temperature
    #fix 3 all ave/time 1 10000 10000 v_v v_t file AgI.vt

    # ---------- Specify ensemble  ---------------------
    fix 4 all npt temp $TEMPERATURE $TEMPERATURE $TDAMP tri 0.0 0.0 1.0

    # --------- Run -------------
    run $NSTEPS
    unfix 4
"""

def compute_AgI_dynamics(timestep, nsteps, temperature, ncpu):
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
                                  intemplate=intemplate, inparam=inparam, ncpu=ncpu)
    output = parse_lammps_thermo(outfile=outfile)
    rdffile = get_rdf(runpath=runpath)
    rdfs = parse_lammps_rdf(rdffile=rdffile)
    return output, rdfs


def md_run():
    output, rdfs = compute_AgI_dynamics(timestep=0.001, nsteps=1000, temperature=300, ncpu=16)
    [simtime, pe, ke, energy, temp, press, dens, msd] = output
    ## ------- plot output properties
    #plt.plot(simtime, temp)
    #plt.show()
    plt.plot(simtime, press)
    plt.show()

    # ----- plot radial distribution functions
    for rdf in rdfs:
        plt.plot(rdf[0], rdf[1])
    plt.show()


if __name__ == '__main__':
    # put here the function that you actually want to run
    md_run()
