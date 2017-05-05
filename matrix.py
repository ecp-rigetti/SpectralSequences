from sage.all import *
from demo import *


class Matrix:
  def __init__(self, R, m, A):
    self.base_ring = R
    self.matrix = matrix(m)
    self.algebra = A

  # Polynomial matrix kernels, brought to you by the magic of syzygies
  def syzygy_kernel_test(self):
    lm = list(self.matrix)
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

  def kernel(self):
    gens = vector(gen_combs(self.base_ring))
    columns = self.syzygy_kernel_test().columns()
    return map(lambda c: gens.dot_product(c), columns)

  def test_boundaries(self):
    lm = list(self.matrix)
    lm = map(list, lm)
    f = lambda lst: str(lst).replace('[', '{').replace(']','}')
    lm = f(lm)
    macaulay2.set('R', 'GF(2)' + str(list(self.base_ring.gens())))
    macaulay2.set('i', 'gens trim image matrix (R,' + lm + ')')
    i = macaulay2('i').to_sage()
    gens = vector(gen_combs(self.base_ring))
    columns = i.columns()
    return map(lambda c: gens.dot_product(c), columns)

  def test_homology(self):
    cycles = self.kernel()
    bounds = self.test_boundaries()
    I = self.base_ring.ideal(cycles)
    reduced = map(I.reduce, bounds)
    print reduced

  # Computational Commutative Algebra, 71
  # Thm 1.6.4
  # monomial_reduce gives us LM(v)/LM(g_i), g_i
  def polynomial_reduction(self, elt):
    R = self.base_ring
    G = self.kernel()
    s = len(G)
    q = [0] * s
    p = 0
    v = elt
    while True:
      while True:
        lm_v = v.lm()
        for g in G:
          g_lm = g.lm()
          try:
            if R.monomial_divides(g_lm, lm_v):
              r = R.monomial_quotient(lm_v, g_lm)
              l = None
              try:
                l = r.exponents()[0]
              except IndexError:
                continue
              degs = list(l)
              evens = map(lambda x: x % 2 == 0, degs)
              alleven = reduce(lambda a, b: a and b, evens)
              if alleven:
                i = G.index(g)
                q[i] += r
                v -= (r * g)
                break
          except ZeroDivisionError:
            continue
        else:
          break
      p += v.lm()
      v -= v.lm()
      if v == 0:
        break
    return (q, p)

  def all_reduce(self):
    bounds = self.test_boundaries()
    return map(self.polynomial_reduction, bounds)









def main():
  # R = PolynomialRing(GF(2), 'x,y,z')
  # R.inject_variables()
  # m = [[x**2,y**2,0],[1,z**2,y**2],[z**2,0,0]]
  # M = Matrix(R, m)
  # print M.matrix
  # print "\n"
  # print M.syzygy_kernel_test()
  h10 = Variable('h10', 1, 0)
  h11 = Variable('h11', 1, 1)
  h12 = Variable('h12', 1, 2)
  h20 = Variable('h20', 2, 0)
  h21 = Variable('h21', 2, 1)
  h30 = Variable('h30', 3, 0)
  gens = [h10,h11,h12,h20,h21,h30]
  A = define_algebra(2, gens, (1,1,1,1,1,1))
  d = first_diff(gens)
  d = parse_diffs(A,d)
  B = define_cdga(A,d)
  m = ap_differential_matrix(B)
  # print_ap_differential_matrix(B)
  R = polynomial_ring_of(A)
  R.inject_variables()
  m = map(lambda x: map(lambda y: eval(special_str(y)), x), m)
  M = Matrix(R,m,B)
  m = M.matrix

  # a = M.syzygy_kernel_test()
  # hopefully_zero = m * a
  # print hopefully_zero.str()
  print M.all_reduce()



if __name__ == "__main__":
  main()
