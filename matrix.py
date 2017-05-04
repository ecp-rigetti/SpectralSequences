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
    print I.reduce(bounds[0])





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
  M.test_homology()



if __name__ == "__main__":
  main()
