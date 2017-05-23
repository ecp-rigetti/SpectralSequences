import sys
import os
from sage.all import *
from functools import *
import itertools


class SpectralSequence:

    """
    Initializes a Spectral Sequence over finite field (of order 2) given
    the base polynomial ring over which the generators are defined, those
    generating elements, the degrees of the generating elements, and a
    dictionary defining the page's differential

    INPUT:

    - ``ring`` -- base polynomial ring (may be a quotient ring) defining the
                                generators of the E1-page of the spectral sequence
    - ``gens`` -- tuple of string representation of generators of this page
    - ``degrees`` -- tuple of degrees of the generators
    - ``diff`` -- dictionary of this page's differential's action on the
                                generators of the page (in string representation)

    """
    def __init__(self, ring, gens, degrees, diff):
        self.page = 1
        self.char = 2
        A = self.define_algebra(self.char, gens, degrees)
        self.covering_ring = self.polynomial_ring_of(A)
        self.relations = []
        self.base_quotient_ring = self.covering_ring.quotient(self.relations)
        self.cover = self.base_quotient_ring.cover()
        self.lift = self.base_quotient_ring.lift()
        self.differential = self.parse_diffs(A, diff)
        self.CDGA = define_cdga(A, self.differential)
        self.generators = self.CDGA.gens()
        self.gencombs = self.gen_combs(A)
        print "Placeholder for the page info declaration"
        print "i.e., planning on listing the pagenum, generators, relations"

    """
    Defines a Graded Commutative Algebra over a Galois Field given characteristic,
    generators, and degrees.

    INPUT:

    - ``char`` -- prime characteristic of Galois Field
    - ``gens`` -- generators of the algebra as a list of objects of type Variable
    - ``degs`` -- tuple of degrees of each generator

    EXAMPLE:

            h10 = Variable('h10', 1, 0)
            h11 = Variable('h11', 1, 1)
            h20 = Variable('h20', 2, 0)
            define_algebra(2, [h10, h11, h20], (1, 1, 2))

    """
    def define_algebra(self, gens, degs):
        g = ''
        for gen in gens:
            g += gen.name + ','
        g = g[:-1]
        A = GradedCommutativeAlgebra(GF(self.char), g, degrees=degs)
        return A

    # Does a function without printing to stdout
    def silence_task(self, f):
        sys.stdout = open(os.devnull, "w")
        f()
        sys.stdout = sys.__stdout__

    """
    Parses a dictionary of differentials represented as strings into sage
    expressions on the algebra A.

    INPUT:

    - ``A`` -- the base algebra, which must contain the variables used in the
                                             differentials
    - ``diffs`` -- a dictionary containing the differentials with string
                                                             representation, such as the output of first_diff

    """
    def parse_diffs(self, A, diffs):
        algebra_diffs = {}
        self.silence_task(A.inject_variables)
        for gen in diffs.keys():
            d = diffs[gen]
            if d == '0':
                    continue
            key = eval(gen)
            value = eval(d)
            algebra_diffs[key] = value
        return algebra_diffs

    """
    Defines a CDGA given a GCA and a dictionary of differentials

    INPUT:

    - ``A`` -- a graded commutative algebra
    - ``diffs`` -- a dictionary of differentials

    EXAMPLE:

            h10 = Variable('h10', 1, 0)
            h11 = Variable('h11', 1, 1)
            h20 = Variable('h20', 2, 0)
            define_algebra(2, [h10, h11, h20], (1, 1, 2))

    """
    def define_cdga(self, A, diffs):
        return A.cdg_algebra(diffs)


    """
    Returns all the sub-A2 combinations of generators of A as a list

    INPUT:

    - ``A`` -- an algebra

    EXAMPLE:

            If A is an algebra with generators h10, h11, h20:

            gen_combs(A) -> [1, h10, h11, h20, h10*h11, h10*h20, h11*h20, h10*h11*h20]

    """
    def gen_combs(self, A):
        t = self.generators
        p = [(1,)]
        l = len(t)
        for i in range(1, l + 1):
            p += itertools.combinations(t, i)
        def mult(tup):
            if len(tup) > 1:
                return tup[0] * mult(tup[1:])
            else:
                return tup[0]
        return map(mult, p)


    # Okay, I guess this is important right now
    # This may be important, but not here
    # IMPORTANT: always inject variables of the corresponding algebra after
    # working with the associated polynomial ring
    def polynomial_ring_of(self, A):
        t = self.generators
        gens = ''
        for gen in t:
            gens += str(gen) + ','
        gens = gens[:-1]
        return PolynomialRing(GF(self.char), gens)


    def special_str(self, x):
        return str(x).replace('^', '**')

    """
    Parses an element of an Algebra as a member of an associated A^2-module with
    prime characteristic p, returning a dictionary of corresponding Ap elements and
    their coefficients which sum to the original element.

    Note: this is a bit of a hack. I couldn't find a gcd function for graded
    commutative algebras, so I managed to do this by converting the elements in the
    algebra to string format, defining a polynomial ring with the same generators as
    the algebra, evaluating the string representations, performing the algorithm,
    converting back to string format, injecting the algebra's generators as
    variables, and finally evaluating the Ap member and the coefficient term as
    members of the algebra.

    INPUT:

    - ``A`` -- an algebra
    - ``p`` -- prime characteristic of the Ap module (A is generally an algebra over
                         GF(p), so this p should be the same)
    - ``term`` -- term in A to be interpreted as product of an Ap term and a
                                coefficient

    EXAMPLE: (outdated, from when this was wrong for sums and returned a tuple)

        Where A =
        Commutative Differential Graded Algebra with generators ('h10', 'h11', 'h20')
        in degrees (1, 1, 1) over Finite Field of size 2 with differential:
         h10 --> 0
         h11 --> 0
         h20 --> h10*h11

        parse_as_Ap(A, 2, boundaries(A).values()[0]) -> (1, h10^2*h11^2)

    """
    def parse_as_A2(self, A, term):
        p = 2
        res = {}
        t = A.gens()
        gens = map(lambda c: c**p, t)
        gens = map(self.special_str, gens)
        term = self.special_str(term)
        R = self.base_quotient_ring

        one = R.one()
        terms = term.split('+')
        terms = map(eval, terms)
        for i in range(0, len(terms)):
            aps = []
            for gen in gens:
                var = eval(gen)
                ap = one
                while True:
                    g = gcd(var, terms[i])
                    if g == var:
                        ap *= g
                        terms[i] /= g
                    else:
                        break
                aps.append(ap)
            ap = one
            for a in aps:
                ap *= a
            if terms[i] in res.keys():
                res[terms[i]] += ap
            else:
                res[terms[i]] = ap
        strres = {}
        for key in res.keys():
            strres[str(key).replace('^', '**')] = str(res[key]).replace('^', '**')
        A.inject_variables()
        final_res = {}
        for key in strres.keys():
            final_res[eval(key)] = eval(strres[key])
        return final_res

    """
    Returns a dictionary of the results of an algebra's differential(s) on the
    coefficient terms of the Ap-module representation of the algebra.

    INPUT:

    - ``A`` -- a CDGA

    EXAMPLE:

        Where A =
        Commutative Differential Graded Algebra with generators ('h10', 'h11', 'h20')
        in degrees (1, 1, 1) over Finite Field of size 2 with differential:
         h10 --> 0
         h11 --> 0
         h20 --> h10*h11

        ap_differential(A) ->
        {h10: (0, 0), 1: (0, 0), h10*h20: (h11, h10^2), h10*h11*h20: (1, h10^2*h11^2),
         h20: (h10*h11, 1), h11*h20: (h10, h11^2), h10*h11: (0, 0), h11: (0, 0)}

    """
    def ap_differential(self, A):
        lst = self.gen_combs(A)
        res = {}
        t = A.gens()
        d = A.differential()
        p = int(str(A.base_ring())[-1:])
        for var in lst:
            dv = d(var)
            if dv == 0:
                res[var] = {0: 0}
            else:
                res[var] = parse_as_Ap(A, p, dv)
        return res


    """
    Returns Sage Matrix representation of of the Ap-module differential, with
    elements considered elements over the Polynomial Ring self.base_ring.

    INPUT:

    - ``A`` -- a CDGA

    EXAMPLE:

        Where A =
        Commutative Differential Graded Algebra with generators ('x', 'y') in degrees
        (1, 2) over Finite Field of size 2 with differential:
         x --> y
         y --> 0

        ap_differential_matrix(A) ->
        [[0, 0, 0, y^2],
         [0, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 0]]

        print_ap_diff_matrix(A) ->
                 1:   x:   y:   x*y:
        1:   0    0    0    y^2
        x:   0    0    0    0
        y:   0    1    0    0
        x*y: 0    0    0    0

    """
    def ap_differential_matrix(self, A):
        lst = gen_combs(A)
        length = len(lst)
        d = ap_differential(A)
        outputs = []
        for i in range(0, length):
            outputs.append([0] * length)
        for i in range(0, length):
            for j in range(0, length):
                try:
                    outputs[i][j] = d[lst[j]][lst[i]]
                except KeyError:
                    outputs[i][j] = 0
        self.silence_task(self.base_quotient_ring.inject_variables)
        m = map(lambda x: map(lambda y: eval(special_str(y)), x), outputs)
        return matrix(m)

    # Polynomial matrix kernels, brought to you by the magic of syzygies!
    def polynomial_matrix_kernel(self, A):
        mat = self.ap_differential_matrix(A):
        lm = list(mat)
        lm = map(lambda x: Ideal(list(x)).syzygy_module().transpose(), lm)
        lm = map(list, lm)
        lm = map(lambda x: map(list, x), lm)
        f = lambda lst: str(lst).replace('[', '{').replace(']','}')
        lm = map(f, lm)
        macaulay2.set('R', 'GF(2)' + str(list(self.base_ring.gens())))
        modules = []
        for i in range(0, len(lm)):
            modules.append('x%d' % i)
            macaulay2.set('x%d' % i, 'image matrix (R,' + lm[i] + ')')
        macaulay2.get('x1')
        modules = str(modules).replace('[','{').replace(']','}').replace('\'', '')
        macaulay2.set('I', 'gens intersect ' + modules)
        # I = macaulay2.get('I')
        I = macaulay2('I').to_sage()
        # print I.str()
        return I


    def turn_page(self, gens, relts, diffs):
            # Make sure to filter differential dictionary to only nonzero diffs
            self.page += 1
            if self.diffs is not {}:








def main:
        pass


if __name__ == "__main__":
        main()
