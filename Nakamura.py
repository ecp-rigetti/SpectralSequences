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

zero = Variable(0, 0)

# Defines a monomial of variables as a list of variables in the monomial
class Monomial:
  def __init__(self, vrs):
    self.vars = vrs
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
    self.s = monomials[0].s
    self.t = monomials[0].t
    self.u = monomials[0].u

  def str_rep(self):
    ret = self.monomials[0].str_rep()
    for i in range(1, len(self.monomials)):
      ret += ' + ' + self.monomials[i].str_rep()
    return ret

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
def first_diff(gens):
  diffs = {}
  for gen in gens:
    i = gen.i
    j = gen.j
    diff = [Monomial([zero])]
    for k in range(1, i):
      if k == 1:
        diff = []
      d1 = next((x for x in gens if x.i == k and x.j == j))
      d2 = next((x for x in gens if x.i == (i - k) and x.j == (k + j)))
      diff.append(Monomial([d1, d2]))
    diffs[gen] = Polynomial(diff)
  return diffs


def main():
  h10 = Variable(1, 0)
  h11 = Variable(1, 1)
  h20 = Variable(2, 0)
  gens = [h10, h11, h20]
  print first_diff(gens)[h20].str_rep()









if __name__ == "__main__":
  main()
