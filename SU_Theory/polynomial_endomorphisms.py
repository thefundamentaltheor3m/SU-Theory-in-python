from numpy import array, Infinity
from sympy import poly, Poly, expand, symbols, sqrt, latex # NOQA F401 [`sqrt` imported for ease of testing]
from numbers import Number


"""Throughout this file, we will assume that our base field has char 0.
Recall that any field of char 0 contains ℚ as a subfield. Hence,
it makes sense to scalar-multiply by any rational number.
"""

x, y, z, a, b, c = symbols('x y z a b c')
_vars = [x, y, z]
_syms = [a, b, c]
# TODO: Implement lexicographic ordering on ℚⁿ.
e1 = array([1, 0, 0])
e2 = array([0, 1, 0])
e3 = array([0, 0, 1])


class DimensionError(ValueError):
    pass


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
            raise DimensionError("Dimension Mismatch: too many polynomials!")
        return poly(expand(s))

    def __mul__(self, G):
        try:
            return Polynomial_Endomorphism(*[
                self(g) for g in G.polys
            ])
        except DimensionError:
            raise DimensionError("Dimension Mismatch: too many polynomials!")

    def __rmul__(self, G):
        if isinstance(G, Number):
            return Polynomial_Endomorphism(*[
                G * p for p in self.polys
            ])
        else:
            return G * self

    def degs(self, weights):
        return array([
            w_degree(p, weights) for p in self.polys
        ])

    def deg(self, weights):
        return sum(self.degs(weights))

    def w_terms(self, weights):
        return array([
            highest_degree_terms(p, weights) for p in self.polys
        ])

    def __getitem__(self, index):
        return self.polys[index].as_expr()

    def __len__(self):
        return len(self.polys)

    def __add__(self, G):
        try:
            return Polynomial_Endomorphism(*[
                self.polys[i] + G.polys[i] for i in range(len(self))
            ])
        except IndexError:
            raise DimensionError("Dimension Mismatch: too many polynomials!")

    def __radd__(self, G):
        return self + G

    def __sub__(self, G):
        return self + (-1 * G)

    def __rsub__(self, G):
        return (-1 * self) + G

    def __iter__(self):
        return self.polys.__iter__()

    def tolatex(self):
        op = "\\begin{bmatrix} " + latex(self[0])
        for i in range(1, len(self)):
            op += " \\\\ " + latex(self[i])
        op += " \\end{bmatrix}"
        return op


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
    if len(w[0]) != len(vars):
        raise ValueError("No. of weights must equal no. of variables!")
    deg = -Infinity
    deg_Vec = [0 for j in range(0,len(vars))]
    terms = _getterms(p)
    for t in terms:
        thisdeg = 0
        thisdegVec = [ 0 for j in range(0,len(vars)) ]
        for i in range(len(vars)):
            thisdeg += w[0][i] * w[1][i] * t[0][i]
            thisdegVec[i] = w[1][i] * t[0][i]
        if deg < thisdeg:
            deg = thisdeg
            deg_Vec = thisdegVec
    return (deg,deg_Vec)


def highest_degree_terms(p: Poly, w: array, vars=_vars):
    res = poly(0, gens=vars)
    terms = _getterms(p)
    d = w_degree(p, w)[0]
    for t in terms:
        this_term = poly(t[1], gens=vars)
        for i in range(len(vars)):
            this_term *= vars[i]**t[0][i]
        if w_degree(this_term, w)[0] == d:
            res += this_term
    return res


# This is for eliminating those terms that are impossible to have the highest degree.

def onlyHigherTerms(P: Poly):
    terms = _getterms(P)
    higher = [i for i in range(0,len(terms))]
    compare = 0
    for i in range(1,len(terms)):
        for j in range(0,i):
            compare = 0
            if terms[i][0][0] < terms[j][0][0]:
                compare = -1
            else:
                compare = 1
            for c in range(1,len(terms[0][0])):
                if (terms[i][0][c] > terms[j][0][c]) and (compare == -1):
                    compare = 0
                    break
                if (terms[i][0][c] < terms[j][0][c]) and (compare == 1):
                    compare = 0
                    break
            if compare == 1:
                higher.remove(j)
                break
            else:
                if compare == -1:
                    higher.remove(i)
                    break
            compare = 0
    theRemaining = []
    for h in higher:
        theRemaining.append(terms[h])
    return theRemaining

def _getterms(p):
    if isinstance(p, Poly):
        terms = p.terms()
    else:
        terms = poly(p).terms()
    return terms
