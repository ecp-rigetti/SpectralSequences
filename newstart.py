from sage.all import *
from functools import *
from moduleops import *
import itertools

"""
with some components inspired by commutative_dga.py from Sage
"""
class CDGA_char2:
    """
    generators writen as list of strings, degrees as tuple
    differential as string dict
    """
    def __init__(self, gens, diff, degrees=None):
        #specify ring definition
        self.p = 2
        self.R = PolynomialRing(GF(2), gens)
        self.MC = Module_Computations(self.R)
        if degrees == None:
            self.degrees = tuple([1] * len(gens))
        else:
            self.degrees = degrees
        #define the differential
        self.R.inject_variables()
        differential = {}
        for gen in self.R.gens():
            differential[gen] = self.R.zero()
        for key in diff.keys():
            differential[self.R(key)] = self.R(diff[key])
        for key in differential.keys():
            val = differential[key]
            if (not val.is_zero()
                    and (not self.is_homogeneous(val) or
                         self.degree(val) != self.degree(key) + 1)):
                raise ValueError("The given dictionary does not determine a degree 1 map")
        self._dic_ = differential
        for gen in self.R.gens():
            if not self.diff(self.diff(gen)).is_zero():
                raise ValueError("The given dictionary does not determine a valid differential")

    #returns degree of an element
    def degree(self, elt):
        if elt.is_zero():
            raise ValueError("The zero element does not have a well-defined degree")
        exps = elt.dict().keys()
        degs = self.degrees
        n = self.R.ngens()
        l = [sum(e[i] * degs[i] for i in range(n)) for e in exps]
        return max(l)

    #checks if well-defined/homogeneous
    def is_homogeneous(self, elt):
        degree = None
        for m in elt.monomials():
            if degree == None:
                degree = self.degree(m)
            else:
                if degree != self.degree(m):
                    return False
        return True

    #calls the differential
    def diff(self, x):
        if x.is_zero():
            return self.R.zero()
        res = self.R.zero()
        dic = x.dict()
        for key in dic:
            keyl = list(key)
            coef = dic[key]
            idx = 0
            while keyl:
                exp = keyl.pop(0)
                if exp > 0:
                    v1 = (exp * self._dic_[x.parent().gen(idx)]
                          * x.parent().gen(idx)**(exp-1))
                    v2 = prod(x.parent().gen(i+idx+1)**keyl[i] for i in
                              range(len(keyl)))
                    res += coef*v1*v2 #semi-neglected coefficient?
                    coef *= x.parent().gen(idx)**exp
                idx += 1
        return res

    # all k-combinations of generators
    def gen_combs(self):
        t = self.R.gens()
        p = [(self.R.one(),)]
        l = len(t)
        for i in range(1, l+1):
            p += itertools.combinations(t, i)
        return map(lambda lst: reduce(lambda x, y: x*y, lst), p)

    # parses algebra element into element of associated A^2 module as a dict
    # with format = basis element (aka gen comb elt): a^2 coefficient
    def parse_as_A2(self, term):
        res = {}
        gens = map(lambda g: g**self.p, self.R.gens())
        if term == self.R.zero():
            return {self.R.zero(): self.R.zero()}
        terms = term.monomials()
        one = self.R.one()
        for i in range(0, len(terms)):
            if terms[i] == 0:
                continue
            aps = []
            for gen in gens:
                ap = one
                while True:
                    g = gcd(gen, terms[i])
                    if g == gen:
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
        return res

    # returns dict of result of algebra's differential on basis elements
    # of ap module
    def ap_differential(self):
        lst = self.gen_combs()
        res = {}
        t = self.R.gens()
        for var in lst:
            dv = self.diff(var)
            res[var] = self.parse_as_A2(dv)
        return res

    #returns a sage matrix rep'n of the ap-module differential
    def ap_differential_matrix(self):
        lst = self.gen_combs() # should really try to remove this, not efficient
        length = len(lst)
        d = self.ap_differential()
        outputs = []
        for i in range(0, length):
            outputs.append([0] * length)
        for i in range(0, length):
            for j in range(0, length):
                try:
                    outputs[i][j] = d[lst[j]][lst[i]]
                except:
                    outputs[i][j] = 0
        return matrix(outputs)

    #computes intersection of two P-submodules (used for intersecting syzygies)
    # takes in two lists of gens
    # following prop 3.2.3 in computational commutative algebra
    def module_intersection(self, gens1, gens2):
        #print 'gens'
        #print gens1
        #print gens2
        all_gens = gens1 + gens2
        syz = self.R.ideal(all_gens).syzygy_module().rows()
        res = []
        n = len(gens1)
        for j in range(len(syz)):
            gen = sum(map(lambda x, y: x*y, syz[j][:n], gens1))
            res.append(gen)
        #print res
        return res


    # calculates cycles of CDGA differential
    def cycles(self):
        gen_combs = vector(self.gen_combs())
        m = self.ap_differential_matrix().columns()
        N = self.MC.syzygy(m)
        cycles = []
        for col in N.columns():
            cycles.append(gen_combs.dot_product(vector(col)))
        return cycles





def main():
    A = CDGA_char2(gens='x, y', degrees=(1, 2) ,diff={'x':'y'})
    print A.cycles()
    #print A.module_intersection([x1,x2], [x1**2-x2**2,x1*x2*x3,x3**2-x1])


if __name__ == "__main__":
    main()
