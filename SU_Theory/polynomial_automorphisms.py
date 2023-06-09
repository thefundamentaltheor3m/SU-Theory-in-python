import numpy as np
import sympy as sp


"""Throughout this file, we will assume that our base field has char 0.
Recall that any field of char 0 contains ℚ as a subfield. Hence,
it makes sense to scalar-multiply by any rational number.
"""


class Polynomial_Automorphism:
    def __init__(self, *polys):
        self.polys = np.array([sp.poly(p) for p in polys])
