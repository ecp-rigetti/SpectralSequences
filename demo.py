from sage.all import *
from functools import *
from matrix import *
import itertools


#TODO: Lump this all into a single "Page" class


#Gradings taken from Ravenel
class Variable:
    def __init__(self, name, i, j):
        self.name = name
        self.i = i
        self.j = j
        #First grading
        self.s = 1
        #Second grading
        self.t = (2**j) * (2**i - 1)
        #May filtration grading (?)
        self.u = i


"""
Calculates the first differential of a May spectral sequence, given the
generators of the first page of the sequence.

INPUT:

- ``gens`` -- generators of the E1-page of the spectral sequence, as a list
                            of Variable objects

OUTPUT: a dictionary whose keys are the string representations of the generators
                of the E1-page and whose values are the string representations of the
                differentials

EXAMPLE:

    h10 = Variable('h10', 1, 0)
    h11 = Variable('h11', 1, 1)
    h20 = Variable('h20', 2, 0)
    first_diff([h10, h11, h20]) -> {'h10': '0', 'h11': '0', 'h20': 'h10h11'}

"""
def first_diff(gens):
    diffs = {}
    for gen in gens:
        i = gen.i
        j = gen.j
        diff = '0'
        for k in range(1, i):
            diff = '' if k == 1 else diff + ' + '
            d1 = next((x for x in gens if x.i == k and x.j == j))
            d2 = next((x for x in gens if x.i == (i - k) and x.j == (k + j)))
            diff += d1.name + '*' + d2.name
        diffs[gen.name] = diff
    return diffs


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
def define_algebra(char, gens, degs):
    g = ''
    for gen in gens:
        g += gen.name + ','
    g = g[:-1]
    A = GradedCommutativeAlgebra(GF(char), g, degrees=degs)
    return A


"""
Parses a dictionary of differentials represented as strings into sage
expressions on the algebra A.

INPUT:

- ``A`` -- the base algebra, which must contain the variables used in the
                     differentials
- ``diffs`` -- a dictionary containing the differentials with string
                             representation, such as the output of first_diff

"""
def parse_diffs(A, diffs):
    algebra_diffs = {}
    A.inject_variables()
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
def define_cdga(A, diffs):
    return A.cdg_algebra(diffs)


"""
Returns all the sub-A2 combinations of generators of A as a list

INPUT:

- ``A`` -- an algebra

EXAMPLE:

    If A is an algebra with generators h10, h11, h20:

    gen_combs(A) -> [1, h10, h11, h20, h10*h11, h10*h20, h11*h20, h10*h11*h20]

"""
def gen_combs(A):
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


# IMPORTANT: always inject variables of the corresponding algebra after
# working with the associated polynomial ring
def polynomial_ring_of(A):
    t = A.gens()
    gens = ''
    for gen in t:
        gens += str(gen) + ','
    gens = gens[:-1]
    p = str(A.base_ring())[-1:]
    return PolynomialRing(GF(p), gens)


"""
Parses an element of an Algebra as a member of an associated Ap-module with
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
def parse_as_Ap(A, p, term):
    res = {}
    t = A.gens()
    gens = map(lambda c: c**p, t)
    gens = map(str, gens)
    for i in range(0, len(gens)):
        gens[i] = gens[i].replace('^', '**')
    term = str(term)
    R = polynomial_ring_of(A)
    R.inject_variables()
    one = R.one()
    term = term.replace('^', '**')
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
def ap_differential(A):
    lst = gen_combs(A)
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
Returns a nested list representation of of the Ap-module differential.
Print function prints this representation.

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
def ap_differential_matrix(A):
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
    return outputs

"""
Prints the Ap-representation differential matrix, with elements tab-delimited

"""
def print_ap_differential_matrix(A):
    outputs = ap_differential_matrix(A)
    outputs = map(lambda l: map(str, l), outputs)
    lst = gen_combs(A)
    top = [' ']
    longest = 0
    for i in range(0, len(lst)):
        el = str(lst[i])
        l = len(el)
        longest = l if l > longest else longest
        outputs[i].insert(0, el + ':')
        top.append(el + ':')
    outputs.insert(0, top)
    for row in outputs:
        for el in row:
            print el.ljust(longest + 1),
        print ''


"""
Helper function that maps a function f to all the keys/values of a dictionary

"""
def dict_map(f, d):
    keys = d.keys()
    out = {}
    for key in keys:
        out[f(key)] = f(d[key])
    return out

"""
Helper function that pulls the gcd of the ap-terms out of a dictionary
representation of the ap-representation of algebra elements.

INPUT:

- ``R`` -- the polynomial ring representation of the base algebra, for gcd use
- ``d`` -- dictionary with ap-representation of algebra element

"""
def pull_a2_gcd(R, d):
    R.inject_variables()
    d = dict_map(special_str, d)
    d = dict_map(eval, d)
    one = R.one()
    keys = d.keys()
    vals = d.values()
    #!!! This is important because python's eval evaluates 1 to an int, not
    #!!! the one element of the polynomial ring
    vals = [one if x == 1 else x for x in vals]
    length = len(keys)
    if length == 1:
        return (d[keys[0]], {keys[0]: 1})
    elif length > 1:
        g = reduce(lambda x, y: gcd(x,y), vals)
        sub = {k: v/g for k, v in d.items()}
        return (g, sub)
    else:
        raise ValueError('Stuff went wrong')

special_str = lambda x: str(x).replace('^', '**')

"""
Note: this comment was for the deprecated version of boundaries; may not be
accurate with output format.

Calculates the boundaries of the Ap-module representation of a CDGA

INPUT:

- ``R`` -- the polynomial ring representation of the base algebra, for gcd use
- ``d`` -- dictionary with Ap-representation of algebra element

OUTPUT: a list of tuples, where elements in the tuple product to the boundary
                elements, the first element of a tuple is the Ap-term, and the second
                element of the tuple is the coefficient term

EXAMPLE:

    Where A =
    Commutative Differential Graded Algebra with generators ('x', 'y') in degrees
    (1, 2) over Finite Field of size 2 with differential:
     x --> y
     y --> 0

    boundaries(A) -> [(y^2, 1), (1, y)]

"""
def boundaries(A):
    diff = ap_differential(A)
    diff_matrix = ap_differential_matrix(A)
    combs = gen_combs(A)
    length = len(diff_matrix)
    totdiff = []
    for d in diff_matrix:
        totdiff += d
    R = polynomial_ring_of(A)
    M = MatrixSpace(R, length, length)
    totdiff = map(special_str, totdiff)
    R.inject_variables()
    totdiff = map(eval, totdiff)
    m = M(totdiff)
    column_indices = m.pivots()
    ap_gcd_with_Ring = partial(pull_a2_gcd, R)
    image = []
    for index in column_indices:
        index = int(str(index))
        A.inject_variables()
        image.append(ap_gcd_with_Ring(diff[combs[index]]))
    f = lambda d: reduce(lambda a, b: a+b, map(lambda (k, v): k * v, d.items()))
    g = lambda (l, r): (l, f(r))
    image = map(g, image)
    image = map(lambda (x, y): (special_str(x), special_str(y)), image)
    A.inject_variables()
    for i in range(0, len(image)):
        a = []
        for j in range(0, len(image[i])):
            a.append(eval(image[i][j]))
        image[i] = tuple(a)
    return image

"""
Prints a nicely-formatted version of the boundaries.

EXAMPLE: (with prior algebra)

print_boundaries(A) -> y^2*A2{1} \oplus 1*A2{y}

"""
def print_boundaries(A):
    p = str(A.base_ring())[-1:]
    bnds = boundaries(A)
    ap = lambda (o, i): str(o) + '*A' + p + '{' + str(i) + '}'
    bnds = map(ap, bnds)
    dir_sum = lambda a, b: a + u' \u2295 ' + b
    print reduce(dir_sum, bnds)


# Incomplete; intended for cycles function
def diff_as_primes(A):
    R = polynomial_ring_of(A)
    R.inject_variables()
    diff = ap_differential_matrix(A)
    diff = map(lambda x: map(special_str, x), diff)
    diff = map(lambda x: map(eval, x), diff)
    primediff = []
    for lst in diff:
        primediff += lst
    t = A.gens()
    t = map(special_str, t)
    t = map(eval, t)
    gens = map(lambda x: x**2, t)
    d = {}
    P = Primes()
    for i in range(0, len(gens)):
        d[gens[i]] = int(P.unrank(i))
    for i in range(0, len(primediff)):
        ac = 1
        for gen in gens:
            while True:
                if primediff[i] == 0 or primediff[i] == 1:
                    break
                else:
                    g = gcd(gen, primediff[i])
                    if g == gen:
                        primediff[i] /= g
                        ac *= d[gen]
                    else:
                        break
        primediff[i] = ac
    return primediff
### TODO: cycles of the algebra
def cycles(A):
    return None



"""
Cohomology functions run on the assumption that the boundaries are subsets
of the cycles (which we haven't yet defined a function to generate)

"""
def cyclic_cohomology_decomp(A):
    bnds = boundaries(A)
    antibnds = filter(lambda (x,y): x == 1, bnds)
    bnds = filter(lambda (x, y): x != 1, bnds)
    print 'Boundaries:'
    print bnds
    print '\"Anti-Boundaries\"'
    print antibnds


def print_cyclic_cohomology(A):
    p = str(A.base_ring())[-1:]
    bnds = boundaries(A)
    ap = lambda (o, i): 'A' + p + '{' + str(i) + '}/' + str(o)
    bnds = map(ap, bnds)
    dir_sum = lambda a, b: a + u' \u2295 ' + b
    print reduce(dir_sum, bnds)


def main():
    # May Example
    # Generators

    h10 = Variable('h10', 1, 0)
    h11 = Variable('h11', 1, 1)
    h12 = Variable('h12', 1, 2)
    h20 = Variable('h20', 2, 0)
    h21 = Variable('h21', 2, 1)
    h30 = Variable('h30', 3, 0)
    gens = [h10, h11, h12, h20, h21, h30]
    # gens = [h10, h11, h20]
    # Defines an algebra
    A = define_algebra(2, gens, (1,1,1,1,1,1))
    # A = GradedCommutativeAlgebra(GF(2), 'x,y,z', (1,1,2))
    d = first_diff(gens)
    d = parse_diffs(A, d)
    B = define_cdga(A, d)

    print_boundaries(B)




    #boundaries(B)
    #print_ap_differential_matrix(B)
    # print_ap_differential_matrix(B)
    #print B.cohomology_generators(6)

    #cyclic_cohomology_decomp(B)
    #print B.cohomology_generators(5)

    # print ap_differential(B)

    # print_boundaries(B)
    # print ""
    # print_cyclic_cohomology(B)
    # print cyclic_cohomology_decomp(B)
    # print B.cohomology_generators(5)

    # Simple Example
    #A = GradedCommutativeAlgebra(GF(2), 'x,y', degrees=(1,2))
    #A.inject_variables()
    #B = A.cdg_algebra({x:y})
    #print ap_differential_matrix(B)

    """
    # Not-so-simple Example
    A = GradedCommutativeAlgebra(GF(2), 'x,y,z', degrees=(1,2,3))
    A.inject_variables()
    B = A.cdg_algebra({x:y, z:y**2})
    #print ap_differential(B)
    #print_boundaries(B)
    #print_boundaries(B)
    #print_ap_differential_matrix(B)
    print_cyclic_cohomology(B)
    #print B
    #print_ap_differential_matrix(B)
    #print boundaries(B)
    #print B.cohomology_generators(10)
    """



if __name__ == "__main__":
    main()


"""
DEPRECATED BOUNDARIES FUNCTION

def boundaries(A):
    res = []
    diff = ap_differential(A).values()
    for val in diff:
        if val != {0: 0}:
            res.append(val)
    dict_str = partial(dict_map, special_str)
    res = map(dict_str, res)
    R = polynomial_ring_of(A)
    R.inject_variables()
    dict_eval = partial(dict_map, eval)
    res = map(dict_eval, res)
    ap_gcd_with_Ring = partial(pull_a2_gcd, R)
    res = map(ap_gcd_with_Ring, res)
    while True:
        length = len(res)
        for i in range(0, length):
            for j in range(i + 1, length):
                #print res
                print i*length + j
                if res[i][1] == res[j][1]:
                    #print res
                    print 'inside'
                    g = gcd(res[i][0], res[j][0])
                    print g
                    if g == res[i][0]:
                        del res[j]
                    elif g == res[j][0]:
                        del res[i]
                    else:
                        res[i] = (res[i][0]/g, res[i][1])
                        res[j] = (res[j][0]/g, res[j][1])
                    break
            else:
                continue
            break
        else:
            break
    f = lambda d: reduce(lambda a, b: a+b, map(lambda (k, v): k * v, d.items()))
    g = lambda (l, r): (l, f(r))
    res = map(g, res)
    res = map(lambda (x, y): (special_str(x), special_str(y)), res)
    A.inject_variables()
    for i in range(0, len(res)):
        a = []
        for j in range(0, len(res[i])):
            a.append(eval(res[i][j]))
        res[i] = tuple(a)
    return res
"""
