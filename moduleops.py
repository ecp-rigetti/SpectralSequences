from sage.all import *
from functools import *

"""
implementing some of the stuff on module grobner basis/syzygies from
kreuzer, computational commutative algebra 1
also taking inspiration from adams and loustaunau, introduction to grobner bases
division algo uses position over term ordering, with lead term being the lead
term of the sage default of the polynomial in lowest i position
"""
class Module_Computations:

    def __init__(self, ring):
        self.R = ring

    """
    gives the leading term of a vector of a module
    uses lead term of generic ordering for lowest index i
    i.e., we look at which has minimal index i with for first element,
    that being the same, look at leading terms
    """
    def lt(self, v):
        if v.is_zero():
            return self.R.zero()
        for i in range(len(v)):
            if v[i].is_zero():
                continue
            vi = self.R(v[i])
            return (vi.lc(), vi.lm(), i)

    """
    determines ordering comparison of two vectors
    returns true if v1 >= v2, false otherwise
    """
    def vec_cmp(self, v1, v2):
        if len(v1) != len(v2):
            raise Exception('not vectors of the same module!')
        lt1 = self.lt(v1)
        lt2 = self.lt(v2)
        if lt1[2] < lt2[2]:
            return True
        elif lt1[2] > lt2[2]:
            return False
        else:
            return lt1[0]*lt1[1] >= lt2[0]*lt2[1]

    """
    runs vectorized division algorithm
    make sure everything has same dimensions!
    """
    def division(self, gs, m):
        gs = filter(lambda g: not g.is_zero(), gs)
        s = len(gs)
        M = FreeModule(self.R, len(gs[0]))
        Mb = M.basis()
        glts = [self.lt(g) for g in gs]
        zero, Mzero = self.R.zero(), M.zero()
        qs = [zero] * s
        p = Mzero
        v = m
        while not v.is_zero(): #pretty jank but somehow works
            b = False
            ltv = self.lt(v)
            for i in range(s):
                ltg_i = glts[i]
                if ltg_i[1].divides(ltv[1]) and ltg_i[2] == ltv[2]:
                    div = (ltv[0]*ltv[1])/(ltg_i[0]*ltg_i[1])
                    #div.reduce() # no clue why this is needed
                    div = self.R(div)
                    qs[i] += div
                    v -= div*gs[i]
                    break
            else:
                e_i = Mb[ltv[2]]
                lt = ltv[0]*ltv[1]*e_i
                p += lt
                v -= lt
            # print v, qs, p
        return (qs, p)


    """
    runs vectorized buchberger algorithm
    """
    def buchberger(self, gs):
        gs = filter(lambda g: not g.is_zero(), gs)
        s = len(gs)
        lts = [self.lt(g) for g in gs]
        B = []
        for i in range(s):
            for j in range(i+1, s):
                if lts[i][2] == lts[j][2]:
                    B.append((i, j))
        while B:
            (i, j) = B.pop()
            lt_i, lt_j = lts[i], lts[j]
            c1 = (lt_j[1] / (lt_i[0] * lt_i[1].gcd(lt_j[1])))
            c2 = (lt_i[1] / (lt_j[0] * lt_i[1].gcd(lt_j[1])))
            # c1.reduce()
            # c2.reduce()
            c1 = self.R(c1)
            c2 = self.R(c2)
            S_ij = (c1 * gs[i]) - (c2 * gs[j])
            div = self.division(gs, S_ij)[1]
            if div.is_zero():
                continue
            else:
                s += 1
                gs.append(div)
                lts.append(self.lt(div))
                for i in range(s):
                    if lts[i][2] == lts[-1][2]:
                        B.append((i, s-1))
        return gs

    """
    also gives a matrix defining how the grobner basis can be expressed
    in terms of original generators
    """
    def ext_buchberger(self, gs):
        gs = filter(lambda g: not g.is_zero(), gs)
        s = len(gs)
        lts = [self.lt(g) for g in gs]
        M = FreeModule(self.R, s)
        A = M.echelonized_basis_matrix().columns()
        B = []
        for i in range(s):
            for j in range(i+1, s):
                if lts[i][2] == lts[j][2]:
                    B.append((i, j))
        while B:
            (i, j) = B.pop()
            lt_i, lt_j = lts[i], lts[j]
            c1 = (lt_j[1] / (lt_i[0] * lt_i[1].gcd(lt_j[1])))
            c2 = (lt_i[1] / (lt_j[0] * lt_i[1].gcd(lt_j[1])))
            # c1.reduce()
            # c2.reduce()
            c1 = self.R(c1)
            c2 = self.R(c2)
            S_ij = (c1 * gs[i]) - (c2 * gs[j])
            div = self.division(gs, S_ij)
            if div[1].is_zero():
                continue
            else:
                s += 1
                gs.append(div[1])
                lts.append(self.lt(div[1]))
                for k in range(s):
                    if lts[k][2] == lts[-1][2]:
                        B.append((k, s-1))
                newcol = (c1*A[i]) - (c2*A[j])
                for q, i in zip(div[0], range(s-1)):
                    newcol -= q*A[i]
                A.append(newcol)
        return gs, matrix(A).transpose()

    """
    computes the syzygy module for a grobner basis
    """
    def grobner_syzygy(self, gs):
        s = len(gs)
        lts = [self.lt(g) for g in gs]
        M = FreeModule(self.R, s)
        basis = M.basis()
        m = []
        B = []
        for i in range(s):
            for j in range(i+1, s):
                if lts[i][2] == lts[j][2]:
                    B.append((i, j))
        while B:
            (i, j) = B.pop()
            lt_i, lt_j = lts[i], lts[j]
            c1 = (lt_j[1] / (lt_i[0] * lt_i[1].gcd(lt_j[1])))
            c2 = (lt_i[1] / (lt_j[0] * lt_i[1].gcd(lt_j[1])))
            # c1.reduce()
            # c2.reduce()
            c1 = self.R(c1)
            c2 = self.R(c2)
            sig_ij = (c1 * basis[i]) - (c2 * basis[j])
            S_ij = (c1 * gs[i]) - (c2 * gs[j])
            if S_ij.is_zero():
                m.append(sig_ij)
            if not S_ij.is_zero():
                fs, r = self.division(gs, S_ij)
                if not r.is_zero():
                    raise Exception('Clearly something went wrong')
                for k in range(s):
                    sig_ij -= fs[k]*basis[k]
                m.append(sig_ij)
        return matrix(m).transpose()

    """
    computes the syzygy module for any system of generators of a submodule
    """
    def syzygy(self, hs):
        gb, A = self.ext_buchberger(hs)
        A = A.rows()
        t = len(hs)
        s = len(gb)
        for i in range(t):
            if hs[i].is_zero():
                A.insert(i, tuple([0] * s))
        A = matrix(A)
        B = []
        for h in hs:
            bs, r = self.division(gb, h)
            B.append(bs)
        B = matrix(B).transpose()
        M = self.grobner_syzygy(gb)
        T = FreeModule(self.R, t)
        I_t = T.echelonized_basis_matrix()
        B_1 = (A*M).columns() if M.ncols() != 0 else [] # M might be empty
        B_2 = (I_t - A*B).columns()
        N = matrix(B_1 + B_2).transpose()
        return N




def main():
    R = PolynomialRing(QQ, 'x,y')
    R.inject_variables()
    M = FreeModule(R, 3)
    gs = [vector([1 - y, 1]), vector([x - y*x, x]), vector([0, 0])]
    # m = vector([y**2 + 2*x**2*y, y**2])
    MC = Module_Computations(R)
    syz = MC.syzygy(gs)
    print syz




if __name__ == "__main__":
    main()

