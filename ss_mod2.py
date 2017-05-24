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
    def __init__(self, gens, degrees, diff):
        self.page = 1
        self.char = 2
        A = self.define_algebra(gens, degrees)
        self.covering_ring = self.polynomial_ring_of(A)
        self.relations = []
        self.base_quotient_ring = self.covering_ring.quotient(self.relations)
        self.cover = lambda x: x
        self.lift = lambda x: x
        self.differential = self.parse_diffs(A, diff)
        self.CDGA = self.define_cdga(A, self.differential)
        self.gencombs = self.gen_combs(A)
        # Just a simple placeholder for var names
        self.n = 0
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
            g += gen + ','
        g = g[:-1]
        A = GradedCommutativeAlgebra(GF(self.char), g, degrees=degs)
        return A

    # Does a function without printing to stdout
    # Very useful for the inject_variables() function
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
        t = A.gens()
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
        t = A.gens()
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
        self.silence_task(R.inject_variables)
        one = R.one()
        terms = term.split('+')
        terms = map(eval, terms)
        terms = map(self.cover, terms)
        for i in range(0, len(terms)):
            aps = []
            for gen in gens:
                var = self.cover(eval(gen))
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
            strres[self.special_str(self.lift(key))] = self.special_str(self.lift(res[key]))
        self.silence_task(A.inject_variables)
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
        for var in lst:
            dv = d(var)
            if dv == 0:
                res[var] = {0: 0}
            else:
                res[var] = self.parse_as_A2(A, dv)
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
        lst = self.gen_combs(A)
        length = len(lst)
        d = self.ap_differential(A)
        outputs = []
        for i in range(0, length):
            outputs.append([0] * length)
        for i in range(0, length):
            for j in range(0, length):
                try:
                    outputs[i][j] = d[lst[j]][lst[i]]
                except KeyError:
                    outputs[i][j] = 0
        self.silence_task(self.covering_ring.inject_variables)
        self.silence_task(self.base_quotient_ring.inject_variables)
        m = map(lambda x: map(lambda y: self.lift(self.cover(eval(self.special_str(y)))), x), outputs)
        return matrix(m)

    # Polynomial matrix kernels, brought to you by the magic of syzygies!
    def kernel(self, m):
        lm = list(m)
        lm = map(lambda x: Ideal(list(x)).syzygy_module().transpose(), lm)
        lm = map(list, lm)
        lm = map(lambda x: map(list, x), lm)
        f = lambda lst: str(lst).replace('[', '{').replace(']','}')
        lm = map(f, lm)
        macaulay2.set('R', 'GF(2)' + str(list(self.covering_ring.gens())))
        macaulay2.set('Q', 'R / ' + f(str(self.relations)))
        modules = []
        for i in range(0, len(lm)):
            modules.append('x%d' % i)
            macaulay2.set('x%d' % i, 'image matrix (Q,' + lm[i] + ')')
        macaulay2.get('x1')
        modules = str(modules).replace('[','{').replace(']','}').replace('\'', '')
        macaulay2.set('I', 'gens intersect ' + modules)
        # I = macaulay2.get('I')
        I = macaulay2('I').to_sage()
        # print I.str()
        gens = vector(self.gen_combs(self.covering_ring))
        columns = I.columns()
        res = map(lambda c: gens.dot_product(c), columns)
        return map(lambda c: self.cover(c), res)

    # Generates boundaries of the differential
    def boundaries(self, m):
        lm = list(m)
        lm = map(list, lm)
        f = lambda lst: str(lst).replace('[', '{').replace(']','}')
        lm = f(lm)
        macaulay2.set('R', 'GF(2)' + str(list(self.covering_ring.gens())))
        macaulay2.set('Q', 'R / ' + f(str(self.relations)))
        macaulay2.set('i', 'gens trim image matrix (Q,' + lm + ')')
        i = macaulay2('i').to_sage()
        gens = vector(self.gen_combs(self.covering_ring))
        columns = i.columns()
        res = map(lambda c: gens.dot_product(c), columns)
        return map(lambda c: self.cover(c), res)

    # Does some weird polynomial reduction on the list of kernel elements mod
    # the kernel and boundary elements to magically produce a multiplicative
    # generating set for the kernel (?)
    # Maybe I should separate (abstract) out the reduction step to clean up code?
    def reduction_step(self):
        ker = self.kernel()[1:]
        b = self.boundaries()
        g = map(lambda x: x**2 , list(self.base_quotient_ring.gens()))
        ker += g
        G = ker + b
        # We begin by addressing the kernel
        while True:
            for i in range(0, len(ker)):
                p = 0
                q = 0
                v = ker[i]
                while True:
                    while True:
                        lm_v = v.lm()
                        r = self.base_quotient_ring.monomial_reduce(lm_v, G[:i] + G[i + 1:])
                        if r == (0, 0):
                            break
                        else:
                            v -= (r[0] * r[1])
                            # It seems that if we can't completely reduce it, we don't want
                            # to reduce it at all (for generating reasons? idk)
                            # So we save this just in case p != 0, to add back in
                            q += (r[0] * r[1])
                    p += v.lm()
                    v -= v.lm()
                    if v == 0:
                        break
                # I don't think this if case actually does anything
                if p != 0:
                    ker[i] = p + q
                    G[i] = p
                else:
                    del ker[i]
                    del G[i]
                    break
            else:
                break
        reduced_kernel = map(lambda c: self.lift(c), ker)
        # Now we address the boundaries
        while True:
            for i in range(0, len(b)):
                p = 0
                q = 0
                v = b[i]
                while True:
                    while True:
                        lm_v = v.lm()
                        r = self.base_quotient_ring.monomial_reduce(lm_v, b[:i] + b[i + 1:])
                        if r == (0, 0):
                            break
                        else:
                            v -= (r[0] * r[1])
                            q += (r[0] * r[1])
                    p += v.lm()
                    v -= v.lm()
                    if v == 0:
                        break
                if p != 0:
                    b[i] = p + q
                else:
                    del b[i]
                    break
            else:
                break
        reduced_boundaries = map(lambda b: self.lift(b), b)
        return (reduced_kernel, reduced_boundaries)

    # This is where the magic happens
    def calculate_generators(self):
        # Make sure to filter differential dictionary to only nonzero diffs
        m = self.ap_differential_matrix(self.CDGA)
        if self.differential is not {}:
            return self.reduction_step(m)
        else:
            return None

    # Gives you a fresh variable name
    def new_varname(self):
        self.n += 1
        return 'var{}'.format(self.n)


    def turn_page(self, gens=None, degs=None, relts=None, diffs=None):
        # This only happens when we have no differentials on a given page
        if gens is None or degs is None or relts is None or diffs is None:
            self.page += 1
        else:
            self.page += 1
            new_gen_names = tuple([self.new_varname() for gen in gens])
            A = self.define_algebra(gens, degs)
            self.covering_ring = self.polynomial_ring_of(A)
            self.relations = []
            self.base_quotient_ring = self.covering_ring.quotient(self.relations)
            self.cover = lambda x: x
            self.lift = lambda x: x
            self.differential = self.parse_diffs(A, diff)
            self.CDGA = self.define_cdga(A, self.differential)
            self.gencombs = self.gen_combs(A)


