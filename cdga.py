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
    def __init__(self, gens, diff, degrees=None, relations=None):
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
        if relations == None:
            self.QR = None
        else:
            relts = [self.R(relt) for relt in relations]
            self.QR = self.R.quotient(relts)
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
        return self.QR.lift(res)

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
        return filter(lambda c: not c.is_zero(), map(lambda c: self.QR.lift(c), cycles))

    # calculates boundaries (image) of CDGA differential
    def boundaries(self):
        gen_combs = vector(self.gen_combs())
        m = self.ap_differential_matrix().columns()
        boundaries = []
        for col in m:
            boundaries.append(gen_combs.dot_product(vector(col)))
        return filter(lambda b: not b.is_zero(), map(lambda b: self.QR.lift(b), boundaries))

    # Does some weird polynomial reduction on the list of kernel elements mod
    # the kernel and boundary elements to magically produce a multiplicative
    # generating set for the cycles and boundaries
    # Maybe I should separate (abstract) out the reduction step to clean up code?
    # Nah.
    # I might have once known how this works; I no longer do.
    def reduction_step(self):
        ker = self.cycles()
        ker.remove(1)
        ker = filter(lambda x: not x.is_zero(), ker)
        b = self.boundaries()
        b = filter(lambda x: not x.is_zero(), b)
        g = map(lambda x: x**2 , self.R.gens())
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
                        v = self.R(v)
                        lm_v = v.lm()
                        r = self.R.monomial_reduce(lm_v, G[:i] + G[i + 1:])
                        print r
                        if r == (0, 0):
                            break
                        elif self.R.monomial_reduce(r[0], ker) == (0, 0): # how the hell does this work???
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
        reduced_kernel = ker # map(lambda c: self.lift(self.cover(c)), ker)
        # Now we address the boundaries
        while True:
            for i in range(0, len(b)):
                p = 0
                q = 0
                v = b[i]
                while True:
                    while True:
                        v = self.R(v)
                        lm_v = v.lm()
                        r = self.R.monomial_reduce(lm_v, b[:i] + b[i + 1:])
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
        reduced_boundaries = b # map(lambda b: self.lift(self.cover(b)), b)
        return (reduced_kernel, reduced_boundaries)



def main():
    A = CDGA_char2(gens='h0, h1, b20', degrees=(1, 1, 2), diff={'b20':'h1**3'}, relations = ['h0*h1'])
    print A.cycles()
    print A.boundaries()
    print A.reduction_step()


if __name__ == "__main__":
    main()
