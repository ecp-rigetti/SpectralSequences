#This will eventually become the main implementation after all the demoing/testing
#is finished in the other implementation. The goal is to combine everything into a
#single spectral sequence page class, which will be able to perpetually progress
#itself through a loop, generating successive pages.

from sage.all import *
from functools import *
import itertools




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


class SS_Page:

  # Initialize the page
  def __init__(self, basis, degrees, char, pagenum):
    self.pagenum = pagenum
    B = self.define_algebra(char, basis, degrees)
    #This assumes a May spectral sequence first page
    string_diff = self.first_diff(basis)
    self.differential = self.parse_diffs(B, string_diff)
    self.A = define_cdga(B, self.differential)
    self.generators = self.A.gens()
    self.gencombs = self.gen_combs(A)


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
  Calculates the first differential of a May spectral sequence, given the
  generators of the first page of the sequence.

  INPUT:

  - ``gens`` -- generators of the E1-page of the spectral sequence, as a list
                of Variable objects

  OUTPUT: a dictionary whose keys are string representations of the generators
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
    self.A.inject_variables()
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


  # IMPORTANT: always inject variables of the corresponding algebra after
  # working with the associated polynomial ring
  def polynomial_ring_of(A):
    t = self.generators
    gens = ''
    for gen in t:
      gens += str(gen) + ','
    gens = gens[:-1]
    p = str(A.base_ring())[-1:]
    return PolynomialRing(GF(p), gens)





def main:
  None


if __name__ == "__main__":
  main()
