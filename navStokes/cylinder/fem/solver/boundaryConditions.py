import numpy as np
from petsc4py import PETSc


class SinuPoiseuilleFlow:

        def __init__(self, t, Um):
            self.t = t
            self.Um = Um

        def set_time(self, t):
            self.t = t

        def __call__(self, x):
            values = np.zeros((2, x.shape[1]), dtype=PETSc.ScalarType)
            values[0] = 4 * self.Um * np.sin(self.t * np.pi / 8) * x[1] * (0.41 - x[1]) / (0.41**2)
            return values

class ConstantPoiseuilleFlow:

        def __init__(self, t, Um):
            self.t = t
            self.Um = Um

        def __call__(self, x):
            values = np.zeros((2, x.shape[1]), dtype=PETSc.ScalarType)
            values[0] = 4 * self.Um * x[1] * (0.41 - x[1]) / (0.41**2)
            return values


