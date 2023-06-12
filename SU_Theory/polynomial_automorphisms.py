from numpy import array
from sympy import poly, Poly, expand, symbols
from sympy.parsing.sympy_parser import parse_expr


"""Throughout this file, we will assume that our base field has char 0.
Recall that any field of char 0 contains ℚ as a subfield. Hence,
it makes sense to scalar-multiply by any rational number.
"""

x, y, z = symbols('x y z')
_vars = [x, y, z]


class Polynomial_Endomorphism:
    def __init__(self, *polys, vars=_vars):
        self.polys = array([poly(expand(p.as_expr()), gens=vars)
                            for p in polys])

    def __call__(self, p: Poly):
        s = str(p.as_expr())
        vars = p.gens
        try:
            for i in range(len(self.polys)):
                s = s.replace(str(vars[i]), f"f{i}")
            for j in range(len(self.polys)):
                s = s.replace(f"f{j}", str(self.polys[j].as_expr()))
        except IndexError:
            raise ValueError("Dimension Mismatch: too many polynomials!")
        return poly(expand(parse_expr(s)))

    def __mul__(self, G):
        try:
            return Polynomial_Endomorphism(*[
                self(g) for g in G.polys
            ])
        except ValueError:
            raise ValueError("Dimension Mismatch: too many polynomials!")


class Polynomial_Automorphism(Polynomial_Endomorphism):
    def __init__(self, *polys):
        super().__init__(self, *polys)
        # TODO: Check invertibility!

    def __call__(self, p: Poly):
        return super().__call__(p)

    def __mul__(self, G):
        if isinstance(G, Polynomial_Automorphism):
            return Polynomial_Automorphism(*((self * G).polys))
        else:
            return super().__mul__(G)
