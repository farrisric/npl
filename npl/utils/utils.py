from copy import deepcopy
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS


def get_crystalline_structure(atoms: Atoms, surrogate_metal: str = 'Cu') -> Atoms:
    original_symbols = deepcopy(atoms.symbols)
    atoms.symbols = [surrogate_metal for _ in atoms.get_positions()]
    atoms.calc = EMT()
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=0.01)
    atoms.symbols = original_symbols
    return atoms
