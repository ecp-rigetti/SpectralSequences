from sage.all import *
from functools import *


# Defines a class for a variable
class Variable:
    def __init__(self, i, j):
        if i != 0 or j != 0:
            self.name = 'h' + str(i) + str(j)
            self.i = i
            self.j = j
            # s, t, u Gradings
            self.s = 1
            self.t = (2**j) * (2**i - 1)
            self.u = i
        else:
            self.name = '0'
            self.i = 0
            self.j = 0
            self.s = 0
            self.t = 0
            self.u = 0

# Defines a monomial of variables as a list of variables in the monomial
class Monomial:
    def __init__(self, vrs):
        self.vars = vrs
        self.name = self.str_rep()
        self.s = sum(v.s for v in self.vars)
        self.t = sum(v.t for v in self.vars)
        self.u = sum(v.u for v in self.vars)

    def str_rep(self):
        ret = self.vars[0].name
        for i in range(1, len(self.vars)):
            ret += '*' + self.vars[i].name
        return ret

# Defines a class for a polynomial of variables
# Note that each monomial must have the same degree
class Polynomial:
    def __init__(self, monomials):
        self.monomials = monomials
        self.name = self.str_rep()
        self.s = monomials[0].s
        self.t = monomials[0].t
        self.u = monomials[0].u

    def is_monomial(self):
        return len(self.monomials) == 1

    def str_rep(self):
        ret = self.monomials[0].str_rep()
        for i in range(1, len(self.monomials)):
            ret += ' + ' + self.monomials[i].str_rep()
        return ret


zero_v = Variable(0, 0)
zero_m = Monomial([zero_v])
zero_p = Polynomial([zero_m])


class Differential:
    # Initialize
    def __init__(self, ring, gens):
        # Base Ring repn. of the CDGA (which I think will remain invariant)
        self.base_ring = ring
        # Store current page
        self.page = 1
        # Store current generators
        self.gens = self.parse_to_vars(gens)
        # Store generators of all pages
        self.all_gens = [tuple(self.gens)]
        # Current differential
        self.differential = self.first_differential()
        # Store differentials of all pages
        self.all_differentials = [self.differential]


    # Parses a list of (i, j) indices as Variables
    def parse_to_vars(self, indices):
        return map(lambda (i, j): Variable(i, j), indices)

    # Parses a list of (i, j) indices as a Monomial
    def parse_to_mon(self, indices):
        return Monomial(self.parse_to_vars(indices))

    # Parses a list of lists of (i, j) indices as a Polynomial
    def parse_to_poly(self, indices_list):
        lst = map(parse_to_mon, indices_list)
        return Polynomial(lst)

    # Parses a list of lists of lists of (i, j) indices as a list of Polynomials
    def parse_to_poly_lst(self, indices_list_list):
        return map(parse_to_poly_lst, indices_list_list)


    # Parses a polynomial as a Sage Polynomial Ring element
    def parse_in_ring(self, poly):
        self.base_ring.inject_variables()
        return eval(poly.str_rep())


    """
    Calculates the first differential of a May spectral sequence, given the
    generators of the first page of the sequence.

    INPUT:

    - ``gens`` -- generators of the E1-page of the spectral sequence, as a list
                                of Variable objects

    OUTPUT: a dictionary whose keys are the variable representation of the
            generators of the E1-page and whose values are the polynomial
            representations of the differentials

    EXAMPLE:

    h10 = Variable(1, 0)
    h11 = Variable(1, 1)
    h20 = Variable(2, 0)
    first_diff([h10, h11, h20]) -> {'h10': '0', 'h11': '0', 'h20': 'h10h11'}

    """
    def first_differential(self):
        diffs = {}
        for gen in self.gens:
            i = gen.i
            j = gen.j
            diff = [zero_m]
            for k in range(1, i):
                if k == 1:
                    diff = []
                d1 = next((x for x in self.gens if x.i == k and x.j == j))
                d2 = next((
                     x for x in self.gens if x.i == (i - k) and x.j == (k + j)))
                diff.append(Monomial([d1, d2]))
            diffs[gen.name] = Polynomial(diff).name
        return diffs


    def kth_diff(self, k):
        diffs = {}
        for gen in self.gens:
            pass


    def turn_page(self, gens):
        self.page += 1
        self.gens = parse_to_poly_lst(gens)
        self.all_gens += [self.gens]
        self.differential = self.kth_diff(self.page)
        self.all_differentials += [self.differential]



def main():
    R = PolynomialRing(GF(2), 'h10,h11,h12,h20,h21,h30')
    R.inject_variables()
    gens = [(1,0),(1,1),(1,2),(2,0),(2,1),(3,0)]
    D = Differential(R, gens)
    print D.differential




if __name__ == "__main__":
    main()
