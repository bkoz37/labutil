from labutil.src.plugins.lammps import *
from ase.spacegroup import crystal
from ase.build import *
import numpy
import matplotlib.pyplot as plt


input_template = """
# ---------- 1. Initialize simulation ---------------------
units metal
atom_style atomic
dimension  3
boundary   p p p
read_data $DATAINPUT

# ---------- 2. Specify interatomic potential ---------------------
pair_style eam/alloy
pair_coeff * * $POTENTIAL  Al

# pair_style lj/cut 4.5
# pair_coeff 1 1 0.392 2.620 4.5

# ---------- 3. Run single point calculation  ---------------------
thermo_style custom step pe lx ly lz press pxx pyy pzz
run 0

# ---- 4. Define and print useful variables -------------
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"

print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
        """

#@wf()
def make_struc(me, alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    unitcell = crystal('Al', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
    #multiplier = numpy.identity(3) * 2
    #ase_supercell = make_supercell(unitcell, multiplier)
    structure = Struc(ase2struc(unitcell))
    return structure


#@wf()
def compute_energy(me, alat, template):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],'Al_zhou.eam.alloy')
    potential = ClassicalPotential(path=potpath, ptype='eam', element=["Al"])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab1", str(alat)))
    struc = make_struc(me, alat=alat)
    output_file = lammps_run(me, struc=struc, runpath=runpath, potential=potential, in_template=template)
    energy, lattice = get_lammps_energy(me, outfile=output_file)
    return energy, lattice


#@wf()
def lattice_scan(me):
    alat_list = numpy.linspace(3.9, 4.3, 6)
    energy_list = [compute_energy(me, alat=a, template=input_template)[0] for a in alat_list]
    print(energy_list)
    plt.plot(alat_list, energy_list)
    plt.show()


if __name__ == '__main__':
    # put here the function that you actually want to run
    #me = Memo(usedb=True)
    b = lattice_scan(me=None)
    #me.finalize()