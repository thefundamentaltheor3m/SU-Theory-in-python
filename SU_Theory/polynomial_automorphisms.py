from numpy import array, Infinity
from sympy import poly, Poly, expand, symbols, sqrt # NOQA F401
import re


"""Throughout this file, we will assume that our base field has char 0.
Recall that any field of char 0 contains ℚ as a subfield. Hence,
it makes sense to scalar-multiply by any rational number.
"""

x, y, z, a, b, c = symbols('x y z a b c')
_vars = [x, y, z]
_syms = [a, b, c]


class Polynomial_Endomorphism:
    def __init__(self, *polys, vars=_vars):
        try:
            polys_as_exprs = [p.as_expr() for p in polys]
        except AttributeError:
            pass
        self.polys = array([poly(expand(p), gens=vars)
                            for p in polys_as_exprs])

    def __call__(self, p: Poly):
        if not isinstance(p, Poly):
            p = poly(p)
        s = p.as_expr()
        vars = p.gens
        try:
            for i in range(len(self.polys)):
                s = s.subs(vars[i], _syms[i])
            for j in range(len(self.polys)):
                s = s.subs(_syms[j], self.polys[j].as_expr())
        except IndexError:
            raise ValueError("Dimension Mismatch: too many polynomials!")
        return poly(expand(s))

    def __mul__(self, G):
        try:
            return Polynomial_Endomorphism(*[
                self(g) for g in G.polys
            ])
        except ValueError:
            raise ValueError("Dimension Mismatch: too many polynomials!")


class Polynomial_Automorphism(Polynomial_Endomorphism):  # Not too important
    def __init__(self, *polys, vars=_vars):
        super().__init__(*polys, vars=vars)
        # TODO: Check invertibility!

    def __call__(self, p: Poly):
        return super().__call__(p)

    def __mul__(self, G):
        # if isinstance(G, Polynomial_Automorphism):
        #     return Polynomial_Automorphism(*((self * G).polys))
        # else:
        #     return super().__mul__(G)
        pass


def w_degree(p: Poly, w: array, vars=_vars):
    if len(w) != len(vars):
        raise ValueError("No. of weights must equal no. of variables!")
    s = ""
    vars = [str(v) for v in vars]
    try:
        s = str(p.as_expr())
    except AttributeError:
        s = str(p)
    deg = -Infinity
    terms = re.split('[+-]', s)
    for t in terms:
        thisdeg = 0
        for c in range(len(t)):
            if t[c] in vars:
                if t[c+1:c+3] == '**':
                    try:
                        thisdeg += int(t[c+3]) * w[vars.index(t[c])]
                    except Exception:
                        raise ValueError("Bad inputs!")
                else:
                    thisdeg += w[vars.index(t[c])]
        deg = max(deg, thisdeg)
    return deg
