"""
Differential class for the mod 2 May spectral sequence

"""
from sage.all import *
from functools import *
from operator import attrgetter, methodcaller


class PolynomialException(Exception):
    pass


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

    def square_0(self):
        self.j += 1

# Defines a monomial of variables as a list of variables in the monomial
class Monomial:
    def __init__(self, vrs):
        if len(vrs) > 1 and zero_v in vrs:
            self.vars = zero_m
        else:
            # We need some sort of sorting to ensure we can effectively
            # decide monomial equality
            self.vars = sorted(vrs, key=attrgetter('i','j'))
        self.name = self.str_rep()
        self.s = sum(v.s for v in self.vars)
        self.t = sum(v.t for v in self.vars)
        self.u = sum(v.u for v in self.vars)

    def is_variable(self):
        l = len(self.vars)
        if l == 0:
            raise PolynomialException("Empty Monomial")
        return l == 1

    # May not actually be a necessary function to have
    def contains_var(self, var):
        return var in self.vars

    # This is to help utilize relations to solve differentials
    def contains_h_1j(self):
        h_1j = filter(lambda x: x.i == 1, self.vars)
        return len(h_1j) != 0

    def str_rep(self):
        ret = self.vars[0].name
        for i in range(1, len(self.vars)):
            ret += '*' + self.vars[i].name
        return ret

# Defines a class for a polynomial of variables
# Note that each monomial must have the same degree
class Polynomial:
    def __init__(self, monomials):
        self.monomials = filter(lambda m: m is not zero_m, monomials)
        if len(self.monomials) == 0:
            self.monomials = [zero_m]
        self.name = self.str_rep()
        self.s = monomials[0].s
        self.t = monomials[0].t
        self.u = monomials[0].u

    def is_monomial(self):
        l = len(self.monomials)
        if l == 0:
            raise PolynomialException("Empty Polynomial")
        return l == 1

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
        self.gens = self.parse_to_poly_lst(gens)
        # Store generators of all pages
        # Seeing how this works as a flattened list for now
        self.all_gens = [tuple(self.gens)]
        # Current differential
        self.differential = self.first_differential()
        # Store differentials of all pages
        self.all_differentials = self.differential
        # Relations (i.e. things quotiented by)
        self.relations = []


    # Parses a list of (i, j) indices as Variables
    def parse_to_vars(self, indices):
        return map(lambda (i, j): Variable(i, j), indices)

    # Parses a list of (i, j) indices as a Monomial
    def parse_to_mon(self, indices):
        return Monomial(self.parse_to_vars(indices))

    # Parses a list of lists of (i, j) indices as a Polynomial
    def parse_to_poly(self, indices_list):
        lst = map(self.parse_to_mon, indices_list)
        return Polynomial(lst)

    # Parses a list of lists of lists of (i, j) indices as a list of Polynomials
    def parse_to_poly_lst(self, indices_list_list):
        return map(self.parse_to_poly, indices_list_list)


    # Parses a polynomial as a Sage Polynomial Ring element
    def parse_in_ring(self, poly):
        self.base_ring.inject_variables()
        return eval(poly.str_rep())

    # Parses a Sage polynomial ring's element into one of my Polynomials
    def parse_from_ring(self, poly):
        p = str(poly)
        p = p.split(' + ')
        p = map(lambda m: m.split('*'), p)
        def expand_mon(mon):
            def expand_expt(s):
                try:
                    i = s.index('^')
                    return [s[:i]] * int(s[i + 1:])
                except ValueError:
                    return [s]
            lst = map(expand_expt, mon)
            return [elt for sub in lst for elt in sub]
        p = map(expand_mon, p)
        def find(elt):
            try:
                return next(var for var in self.all_gens[0] if var.monomials[0].vars[0].name == elt)
            except StopIteration:
                return zero_p
        p = map(lambda lst: map(find, lst), p)
        p = map(lambda lst: map(lambda elt: elt.monomials[0].vars[0], lst), p)
        p = map(lambda lst: Monomial(lst), p)
        print Polynomial(p).name
        return Polynomial(p)

    def poly_prod(self, p1, p2):
        p1, p2 = self.parse_in_ring(p1), self.parse_in_ring(p2)
        p = p1 * p2
        return self.parse_from_ring(p)

    def poly_sum(self, p1, p2):
        p1, p2 = self.parse_in_ring(p1), self.parse_in_ring(p2)
        p = p1 + p2
        return self.parse_from_ring(p)


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
        gens = map(lambda x: x.monomials[0].vars[0], self.gens)
        for genp in self.gens:
            gen = genp.monomials[0].vars[0]
            i = gen.i
            j = gen.j
            diff = [zero_m]
            for k in range(1, i):
                if k == 1:
                    diff = []
                d1 = next((x for x in gens if x.i == k and x.j == j))
                d2 = next((
                     x for x in gens if x.i == (i - k) and x.j == (k + j)))
                diff.append(Monomial([d1, d2]))
            diffs[genp] = Polynomial(diff)
        return diffs

    def square_n(self, n, poly):
        # Checking if the monomial is a single variable via is_variable()
        # guarantees we have a non-empty list, since that function raises
        # an exception at the empty list.
        if poly.is_monomial():
            mon = poly.monomials[0]
            if mon.is_variable():
                if n == 0:
                    print mon.name
                    print mon.vars[0].name
                    mon.vars[0].square_0()
                    t = Polynomial([Monomial([mon.vars[0]])])
                    return t if t in self.all_gens[0] else zero_p
                elif n == 1:
                    t = mon.vars[0]
                    return Polynomial([Monomial([t, t])])
                else:
                    return zero_p
            else:
                x, y = mon.vars[:1], mon.vars[1:]
                x, y = Polynomial([Monomial(x)]), Polynomial([Monomial(y)])
                ret = []
                for i in range(0, n):
                    ret += self.poly_prod(self.square_n(i, x), self.square_n(n - i, y))
                return reduce(self.poly_sum, ret)
        else:
            return reduce(self.poly_sum, map(lambda m: self.square_n(n, m), poly))

    def square_term_differentials(self):
        diffs = {}
        for gen in self.gens:
            if (gen.is_monomial()
                    and gen.monomials[0].is_variable()
                    and gen.monomials[0].vars[0].i == 1):
                continue
            new_rep = self.parse_in_ring(gen)
            squared_gen = new_rep ** 2
            squared_gen = self.parse_from_ring(squared_gen)
            unsquared = self.differential[gen]
            # Completely unfounded assumption that the squaring operator
            # splits over addition
            print unsquared.name
            squared = self.square_n(1, unsquared)
            diffs[squared_gen] = squared
        self.all_differentials.update(diffs)


    # Calculates the result of applying the rth May differential to a
    # Polynomial input; intended to be applied to the generators of a page
    def rth_diff(self, r, poly):
        print poly.name
        if poly.is_monomial():
            mon = poly.monomials[0]
            # If it's just a single variable, it has to be something of form
            # h_1j, because none of the higher terms survive past the first
            # page (since they all have defined differential).
            if mon.is_variable():
                return zero_p
            # If it's a monomial that contains one of the h_1j, we can use
            # a differential trick like this with one of the relations,
            # since kernel elements pull out of the differential:
            # h_1k * d(h_1j * m)         (where k = j +/- 1)
            # = h_1k * h_1j * d(m)
            # = 0                        (which means d(h_1j * m) = 0)
            elif mon.contains_h_1j():
                return zero_p
            # Otherwise, it should be one of the squared terms
            else:
                return self.all_differentials[poly]
        # We otherwise have an unfortunate polynomial to deal with
        else:
            if not any(map(lambda m: m.contains_h_1j(), poly.monomials)):
                raise PolynomialException("You're basically screwed")
            else:
                return zero_p


    def turn_page(self, gens, relations):
        self.page += 1
        self.relations += relations
        self.square_term_differentials()
        self.gens = self.parse_to_poly_lst(gens)
        self.all_gens += [tuple(self.gens)]
        self.differential = self.rth_diff_on_gens(self.page)
        print self.differential
        self.all_differentials += [self.differential]



def main():
    R = PolynomialRing(GF(2), 'h10,h11,h12,h20,h21,h30')
    R.inject_variables()
    gens = [[[(1,0)]],[[(1,1)]],[[(1,2)]],[[(2,0)]],[[(2,1)]],[[(3,0)]]]
    D = Differential(R, gens)
    gens_2 = [h10, h11, h12, h20**2, h21**2, h30**2, h20*h21 + h11*h30]
    relts = [h10*h11, h11*h12, h20*h12 + h21*h10]
    D.turn_page(map(D.parse_from_ring, gens_2), map(D.parse_from_ring, relts))




if __name__ == "__main__":
    main()
