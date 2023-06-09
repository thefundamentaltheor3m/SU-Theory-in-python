from numpy import array
from sympy import poly, Poly, simplify
from sympy.parsing.sympy_parser import parse_expr


"""Throughout this file, we will assume that our base field has char 0.
Recall that any field of char 0 contains ℚ as a subfield. Hence,
it makes sense to scalar-multiply by any rational number.
"""


class Polynomial_Endomorphism:
    def __init__(self, *polys):
        self.polys = array([poly(p) for p in polys])

    def __call__(self, p: Poly):
        s = str(p)
        vars = p.gens
        try:
            for i in range(len(self.polys)):
                s.replace(str(vars[i]), f"({str(self.polys[i])})")
        except IndexError:
            raise ValueError("Dimension Mismatch: too many polynomials!")
        return poly(simplify(parse_expr(s)))


def Jacobian(F):  # The Jacobian Determinant of a polynomial endomorphism
    pass


class Polynomial_Automorphism(Polynomial_Endomorphism):
    def __init__(self, *polys):
        super().__init__(self, *polys)
