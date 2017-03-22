from src.plugins.pwscf import *
from ase.spacegroup import crystal
from ase.io import write
from ase.build import *
import numpy
import matplotlib.pyplot as plt


#@wf()
def make_struc(me, alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    fecell = bulk('Fe', 'bcc', a=alat)
    # check how your cell looks like
    #write('s.cif', gecell)
    structure = Struc(ase2struc(fecell))
    return structure


#@wf()
def compute_energy(me, alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potname = 'Fe.pbe-nd-rrkjus.UPF'
    potpath = os.path.join(os.environ['ESPRESSO_PSEUDO'], potname)
    pseudopots = {'Fe': PseudoPotential(path=potpath, ptype='uspp', element='Fe',
                                        functional='LDA', name=potname)}
    struc = make_struc(me, alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    dirname = 'Fe_a_{}_ecut_{}_nk_{}'.format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Problem1", dirname))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
            'pseudo_dir': os.environ['ESPRESSO_PSEUDO'],
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 8,
            'nspin': 2,
            'starting_magnetization(1)': 0.7,
            'occupations': 'smearing',
            'smearing': 'mp',
            'degauss': 0.02
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
        })

    output_file = run_qe_pwscf(me, runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, ncpu=2)
    output = parse_qe_pwscf_output(me, outfile=output_file)
    return output


#@wf()
def lattice_scan(me):
    nk = 3
    ecut = 30
    alat = 3.0
    output = compute_energy(me, alat=alat, ecut=ecut, nk=nk)
    print(output)

if __name__ == '__main__':
    # put here the function that you actually want to run
    lattice_scan(me=None)
