import ase, json
from .util import *


class Generic(object):
    """
    The basic class for a data object. it is defined by the attribute content which can be a number,
    string, list or a dict.
    """
    def __init__(self, *args, **kwargs):
        if args:
            # in case of a single literal or a dict being supplied
            assert not kwargs and len(args) == 1, "wrong number of init args"
            self.content = args[0]
        if kwargs:
            # can supply content=dict, or an unpacked dict by keywords, but not both
            if 'content' in kwargs:
                assert not args and len(kwargs) == 1, "wrong number of init kwargs"
                self.content = kwargs['content']
            else:
                self.content = kwargs

    def __repr__(self):
        return "{} {}".format(self.__class__.__name__, self.content)


class Path(Generic):
    """data class used to store OS path objects as a dict {'path': <path_string>} """
    @property
    def path(self):
        return self.content['path']


class Dir(Path):
    """data class used to store OS directories"""
    pass


class File(Path):
    """data class used to store OS files"""
    pass


class TextFile(File):

    @property
    def text(self):
        return self.content['text']

    def write(self):
        write_file(fname=self.path, text=self.text)


class ExternalCode(File):
    """Class representing a binary code executable, such as LAMMPS or Quantum ESPRESSO"""
    pass


class Param(Generic):
    """Class representing a dictionary of parameters"""
    pass


class Constraint(Param):
    """
    Data class describing constraints on the structure for relaxation and dynamics calculations
    example:
        constraint = {'atoms': {'3': [0, 0, 0],
                                '4': [1, 1, 1]}
                      'cell' : 'volume'
                      }

    """
    pass


class Struc(Param):
    """
    Data class containing information about a structure
    example:
        struc1 = {"cell": [[1.0, 0, 0],[0, 1.0, 0],[0, 0, 2.0]],
                  "periodicity" : [True, True, True],
                  "species": {'H': {'id': 1, 'mass':1.008, 'atomic_number': 1},
                              'He',{'id': 2, 'mass': 4.003, 'atomic_number' : 2}}
                  "positions": [['H', [4.0, 3.0, 6.0]],
                                ['He', [4.0, 5.0, 9.0]]],
                 }
    """
    @staticmethod
    def from_ase(aseobj):
        # need to use method tolist() of numpy arrays to get valid json
        cell = aseobj.cell.tolist()
        pbc = aseobj.get_pbc().tolist()
        symbols = aseobj.get_chemical_symbols()
        masses = aseobj.get_masses()
        positions = aseobj.get_positions().tolist()
        # easy way to get rid of tuples after zipping
        positions = json.loads(json.dumps(list(zip(symbols, positions))))
        types = sorted(list(set(zip(symbols, masses))))
        species = {tp[0]: {'mass': tp[1], 'kind': i + 1} for i, tp in enumerate(types)}
        content = {'cell': cell, 'positions': positions, 'periodicity': pbc, 'species': species}
        return content

    def to_ase(self):
        cell = self.content['cell']
        atoms = [ase.Atom(site[0], tuple(site[1])) for site in self.content['positions']]
        aseobj = ase.Atoms(atoms)
        aseobj.set_cell(cell)
        aseobj.set_pbc(self.content['periodicity'])
        return aseobj

    @property
    def symbols(self):
        return [s[0] for s in self.content['positions']]

    @property
    def cell(self):
        return self.content['cell']

    @property
    def positions(self):
        return self.content['positions']

    @property
    def n_atoms(self):
        return len(self.content['positions'])

    @property
    def n_species(self):
        return len(self.species)

    @property
    def species(self):
        return self.content['species']


class Kpoints(Param):
    """
    Class describing the Brillouin Zone k-mesh for plane-wave DFT calculations.
    Example:
        kpts = {'option': 'automatic', 'gridsize': [2, 3, 3]}
    """
    pass


class Potential(File):
    """Data class to store information about a Pseudo Potential file"""
    pass


class PseudoPotential(Potential):
    """Data class to store information about a Pseudo Potential file"""
    def __init__(self, **kwargs):
        if 'path' not in kwargs:
            name = kwargs['name']
            potpath = os.path.join(os.environ['ESPRESSO_PSEUDO'], name)
            kwargs.update({'path': potpath})
        super().__init__(**kwargs)


class ClassicalPotential(Potential):
    """ Classical potential, e.g. EAM, Buckingham, OPLS etc """
    def __init__(self, **kwargs):
        if 'path' not in kwargs:
            potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'], kwargs['name'])
            kwargs.update({'path': potpath})
        super().__init__(**kwargs)
    pass


def ase2struc(ase_atoms):
    return Struc.from_ase(ase_atoms)


def struc2ase(struc):
    return Struc.to_ase(struc)


#def wf(func):
#    return func