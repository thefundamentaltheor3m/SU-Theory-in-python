from numpy import array
from sympy import poly, Poly, simplify, symbols
from sympy.parsing.sympy_parser import parse_expr


"""Throughout this file, we will assume that our base field has char 0.
Recall that any field of char 0 contains ℚ as a subfield. Hence,
it makes sense to scalar-multiply by any rational number.
"""

x, y, z = symbols('x y z')
_vars = [x, y, z]


class Polynomial_Endomorphism:
    def __init__(self, *polys, vars=_vars):
        self.polys = array([poly(p, gens=vars) for p in polys])

    def __call__(self, p: Poly):
        s = str(p.as_expr())
        vars = p.gens
        try:
            for i in range(len(self.polys)):
                s = s.replace(str(vars[i]), f"f{i}")
            s = s.replace("f1", str(self.polys[1].as_expr()))
            s = s.replace("f2", str(self.polys[2].as_expr()))
            s = s.replace("f3", str(self.polys[3].as_expr()))
        except IndexError:
            raise ValueError("Dimension Mismatch: too many polynomials!")
        return poly(simplify(parse_expr(s)))


def Jacobian(F):  # The Jacobian Determinant of a polynomial endomorphism
    pass


class Polynomial_Automorphism(Polynomial_Endomorphism):
    def __init__(self, *polys):
        super().__init__(self, *polys)
