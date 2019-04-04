import numpy
from string import Template
from labutil.src.objects import *
from ase.io.lammpsrun import read_lammps_dump


def write_lammps_data(struc, runpath):
    """Make LAMMPS struc data"""
    datatxt = 'Header of the LAMMPS data file \n\n'
    datatxt += '{} atoms \n'.format(struc.n_atoms)
    datatxt += '{} atom types\n\n'.format(struc.n_species)

    datatxt += '0.0 {}  xlo xhi\n'.format(struc.cell[0][0])
    datatxt += '0.0 {}  ylo yhi\n'.format(struc.cell[1][1])
    datatxt += '0.0 {}  zlo zhi\n'.format(struc.cell[2][2])

    datatxt += '\nMasses \n\n'
    for sym, sp in struc.species.items():
        datatxt += '{}  {}\n'.format(sp['kind'], sp['mass'])
    # Write atom positions in angstrom
    datatxt += '\nAtoms # atomic \n\n'
    for index, site in enumerate(struc.positions):
        datatxt += '{} {} {:1.5f} {:1.5f} {:1.5f} \n'.format(index + 1, struc.species[site[0]]['kind'], *site[1])

    datafile = os.path.join(runpath.path, 'lammps.data')
    write_file(datafile, datatxt)
    return File(path=datafile)


def write_lammps_input(datafile, runpath, rdffile, intemplate, inparam, potential):
    """make Lammps input script"""
    if potential:
        ppath = potential.path
    else:
        ppath = ''

    dumpfile = File(path=os.path.join(runpath.path, 'lammps.dump'))
    subst = {'DATAINPUT': datafile.path, 'POTENTIAL': ppath, 'DUMPFILE': dumpfile.path, 'RDFFILE': rdffile.path}
    inptxt = Template(intemplate).safe_substitute({**subst, **inparam})
    infile = TextFile(path=os.path.join(runpath.path, 'lammps.in'), text=inptxt)
    infile.write()
    return infile


def lammps_run(struc, runpath, intemplate, potential, inparam):
    lammps_code = ExternalCode(path=os.environ['LAMMPS_COMMAND'])
    prepare_dir(runpath.path)
    logfile = File(path=os.path.join(runpath.path, 'lammps.log'))
    outfile = File(path=os.path.join(runpath.path, 'lammps.out'))
    rdffile = File(path=os.path.join(runpath.path, 'lammps.rdf'))

    datafile = write_lammps_data(struc=struc, runpath=runpath)
    infile = write_lammps_input(datafile=datafile, potential=potential, runpath=runpath, rdffile=rdffile,
                                intemplate=intemplate, inparam=inparam)

    lammps_command = "{} -in {} -log {} > {}".format(lammps_code.path, infile.path,
                                                     logfile.path, outfile.path)
    run_command(lammps_command)
    return outfile


def get_rdf(runpath):
    rdffile = File(path=os.path.join(runpath.path, 'lammps.rdf'))
    return rdffile


def get_lammps_energy(outfile):
    energy = None
    lattice = None
    with open(outfile.path, 'r') as fout:
        for line in fout.readlines():
            if 'Total energy' in line:
                energy = float(line.split(" = ")[1])
            if 'Lattice constant (Angstoms)' in line:
                lattice = float(line.split(" = ")[1])
    return energy, lattice


def parse_lammps_thermo(outfile):
    """
    Parser for the main LAMMPS output file, extracting columns of the thermo output
    """
    with open(outfile.path, 'r') as fout:
        read_temp = False
        output = []
        for line in fout:
            if not read_temp and line.startswith('Step'):
                read_temp = True
                continue
            if read_temp and line.startswith('Loop'):
                read_temp = False
                continue
            if read_temp:
                values = [float(value) for value in line.split()]
                output.append(values)
    outrows = numpy.transpose(numpy.array(output))
    return outrows


def parse_lammps_rdf(rdffile):
    """
    Parse the RDF file written by LAMMPS
    """
    with open(rdffile.path, 'r') as rdfout:
        rdfs = []
        buffer = []
        for line in rdfout:
            values = line.split()
            if line.startswith('#'):
                continue
            elif len(values) == 2:
                nbins = values[1]
            else:
                buffer.append([float(values[1]), float(values[2])])
                if len(buffer) == int(nbins):
                    frame = numpy.transpose(numpy.array(buffer))
                    rdfs.append(frame)
                    buffer = []
    return rdfs


def parse_structure_dump(runpath, dumpfilename):
    dumpfile = os.path.join(runpath.path, dumpfilename)
    last_structure = read_lammps_dump(dumpfile)
    print(last_structure)
    strucfile = os.path.join(runpath.path, 'struc.cif')
    ase.io.write(strucfile, last_structure)
