from src.lammps import *
from ase.spacegroup import crystal
from ase.build import *
import numpy
import matplotlib.pyplot as plt


input_template = """
# ---------- Initialize simulation ---------------------
units metal
atom_style atomic
dimension  3
boundary   p p p
read_data $DATAINPUT

# ---------- Specify interatomic potential ---------------------
#pair_style eam/alloy
#pair_coeff * * $POTENTIAL  Al

pair_style lj/cut 2.5
pair_coeff 1 1 0.392 2.620 2.5

# ---------- 3. Run the calculation ----------------
# -- perform a single-point energy calculation only
run 0

# -- include optimization of the unit cell parameter
# fix 1 all box/relax iso 0.0 vmax 0.001

# -- enable optimization of atomic positions (and the cell)
# min_style cg
# minimize 1e-10 1e-10 1000 10000

# ---- Define and print useful variables -------------
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"

print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
        """


def make_struc(alat):
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


def compute_energy(alat, template):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potential = ClassicalPotential({'path': os.path.join(os.environ['LAMMPS_POTENTIALS'],'Al_zhou.eam.alloy')})
    runpath = Dir({'path': os.path.join(os.environ['WORKDIR'], "Lab1", str(alat))})
    struc = make_struc(alat=alat)
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, in_template=template)
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice


def lattice_opti():
    a = 3.9
    result = compute_energy(alat=a, template=input_template)
    print(result)


if __name__ == '__main__':
    # put here the function that you actually want to run
    lattice_opti()