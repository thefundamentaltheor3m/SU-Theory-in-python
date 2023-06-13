from numpy import array, Infinity
from sympy import poly, Poly, expand, symbols, sqrt # NOQA F401
from sympy.parsing import parse_expr
import re


"""Throughout this file, we will assume that our base field has char 0.
Recall that any field of char 0 contains ℚ as a subfield. Hence,
it makes sense to scalar-multiply by any rational number.
"""

x, y, z, a, b, c = symbols('x y z a b c')
_vars = [x, y, z]
_syms = [a, b, c]
e1 = array([1,0,0])
e2 = array([0,1,0])
e3 = array([0,0,1])


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

    def deg(self, weights):
        return array([
            w_degree(p, weights) for p in self.polys
        ])


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


def verified_str(p):
    s = ""
    try:
        s = str(p.as_expr())
    except AttributeError:
        s = str(p)
    return s


def w_degree(p: Poly, w: array, vars=_vars):
    if len(w) != len(vars):
        raise ValueError("No. of weights must equal no. of variables!")
    vars = [str(v) for v in vars]
    deg = -Infinity
    s = verified_str(p)
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


def highest_degree_terms(p: Poly, w: array, vars=_vars):
    res = poly(0, gens=vars)
    s = verified_str(p)
    terms = re.split('[+-]', s)
    terms = [parse_expr(t) for t in terms]
    d = w_degree(p, w)
    for t in terms:
        if w_degree(t, w) == d:
            res += t
    return res
