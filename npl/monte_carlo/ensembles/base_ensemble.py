import os
import random
import warnings

from abc import ABC, abstractmethod
from collections import OrderedDict
from math import gcd
from time import time
from typing import Any, Dict, List, Optional, Type, Union

import numpy as np

from ase import Atoms
from ase.calculators.calculator import Calculator

class BaseEnsemble(ABC):
    def __init__(self,
                 structure: Atoms,
                 calculator: Calculator,
                 user_tag: str = None,
                 random_seed: int = None,
                 trajectory_write_interval: int = None) -> None:
        
        self._accepted_trials = 0
        self._step = 0

        # calculator and configuration
        self._calculator = calculator
        self._user_tag = user_tag

        # random number generator
        if random_seed is None:
            self._random_seed = random.randint(0, int(1e16))
        else:
            self._random_seed = random_seed
        random.seed(a=self._random_seed)

    @property
    def structure(self) -> Atoms:
        """ Current configuration (copy). """
        return self.structure
    
    @property
    def step(self) -> int:
        """ Current trial step counter. """
        return self._step




