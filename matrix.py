from sage.all import *
from demo import *

###Dear God have mercy on my soul for attempting this###
class Matrix:
  def __init__(self, R, m):
    self.base_ring = R
    self.matrix = matrix(m)

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
    macaulay2.set('I', 'intersect ' + modules)
    print macaulay2('I_{1}')
    return I








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
  print_ap_differential_matrix(B)
  R = polynomial_ring_of(A)
  R.inject_variables()
  m = map(lambda x: map(lambda y: eval(special_str(y)), x), m)
  M = Matrix(R,m)
  a = M.syzygy_kernel_test()


if __name__ == "__main__":
  main()
