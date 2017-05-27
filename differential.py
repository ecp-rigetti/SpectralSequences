"""
This is a copy of the Differential class in /sage/algebras/commutative_dga.py
which has been modified so that the morphism can be defined on an element which
does not generate the algebra A

"""
class Differential(UniqueRepresentation, Morphism):

    def __init__(self, A, diff):
        r"""
        Initialize ``self``.

        INPUT:

        - ``A`` -- algebra where the differential is defined

        - ``im_gens`` -- differential containing the image of each generator

        EXAMPLES::

            sage: A.<x,y,z,t> = GradedCommutativeAlgebra(QQ)
            sage: B = A.cdg_algebra({x: x*y, y: x*y, z: z*t, t: t*z})
            sage: [B.cohomology(i).dimension() for i in range(6)]
            [1, 2, 1, 0, 0, 0]
            sage: d = B.differential()

        We skip the category test because homsets/morphisms aren't
        proper parents/elements yet::

            sage: TestSuite(d).run(skip="_test_category")

        An error is raised if the differential `d` does not have
        degree 1 or if `d \circ d` is not zero::

            sage: A.<a,b,c> = GradedCommutativeAlgebra(QQ, degrees=(1,2,3))
            sage: A.cdg_algebra({a:b, b:c})
            Traceback (most recent call last):
            ...
            ValueError: The given dictionary does not determine a valid differential
        """
        self._dic_ = diff
        Morphism.__init__(self, Hom(A, A, category=Modules(A.base_ring())))

        for i in A.gens():
            if not self(self(i)).is_zero():
                raise ValueError("The given dictionary does not determine a valid differential")

    def _call_(self, x):
        r"""
        Apply the differential to ``x``.

        INPUT:

        - ``x`` -- an element of the domain of this differential

        EXAMPLES::

            sage: A.<x,y,z,t> = GradedCommutativeAlgebra(QQ)
            sage: B = A.cdg_algebra({x: x*y, y: x*y, z: z*t, t: t*z})
            sage: D = B.differential()
            sage: D(x*t+1/2*t*x*y) # indirect doctest
            -1/2*x*y*z*t + x*y*t + x*z*t

        Test positive characteristic::

            sage: A.<x,y> = GradedCommutativeAlgebra(GF(17), degrees=(2,3))
            sage: B = A.cdg_algebra(differential={x:y})
            sage: B.differential()(x^17)
            0
        """
        if x.is_zero():
            return self.codomain().zero()
        res = self.codomain().zero()
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
                    res += coef*v1*v2
                    coef *= ((-1) ** total_degree(x.parent()._degrees[idx])
                             * x.parent().gen(idx)**exp)
                idx += 1
        return res

    def _repr_defn(self):
        """
        Return a string showing where ``self`` sends each generator.

        EXAMPLES::

            sage: A.<x,y,z,t> = GradedCommutativeAlgebra(QQ)
            sage: B = A.cdg_algebra({x: x*y, y: x*y, z: z*t, t: t*z})
            sage: D = B.differential()
            sage: print(D._repr_defn())
            x --> x*y
            y --> x*y
            z --> z*t
            t --> -z*t
        """
        return '\n'.join("{} --> {}".format(i, self(i)) for i in self.domain().gens())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: A.<x,y,z,t> = GradedCommutativeAlgebra(QQ)
            sage: D = A.differential({x: x*y, y: x*y, z: z*t, t: t*z})
            sage: D
            Differential of Graded Commutative Algebra with generators ('x', 'y', 'z', 't') in degrees (1, 1, 1, 1) over Rational Field
              Defn: x --> x*y
                    y --> x*y
                    z --> z*t
                    t --> -z*t
        """
        if self.domain() is None:
            return "Defunct morphism"

        s = "Differential of {}".format(self.domain()._base_repr())
        s += "\n  Defn: " + '\n        '.join(self._repr_defn().split('\n'))
        return s

    @cached_method
    def differential_matrix(self, n):
        r"""
        The matrix that gives the differential in degree ``n``.

        INPUT:

        - ``n`` -- degree

        EXAMPLES::

            sage: A.<x,y,z,t> = GradedCommutativeAlgebra(GF(5), degrees=(2, 3, 2, 4))
            sage: d = A.differential({t: x*y, x: y, z: y})
            sage: d.differential_matrix(4)
            [0 1]
            [2 0]
            [1 1]
            [0 2]
            sage: A.inject_variables()
            Defining x, y, z, t
            sage: d(t)
            x*y
            sage: d(z^2)
            2*y*z
            sage: d(x*z)
            x*y + y*z
            sage: d(x^2)
            2*x*y
        """
        A = self.domain()
        dom = A.basis(n)
        cod = A.basis(n+1)
        cokeys = [a.lift().dict().keys()[0] for a in cod]
        m = matrix(A.base_ring(), len(dom), len(cod))
        for i in range(len(dom)):
            im = self(dom[i])
            dic = im.lift().dict()
            for j in dic.keys():
                k = cokeys.index(j)
                m[i,k] = dic[j]
        m.set_immutable()
        return m

    def coboundaries(self, n):
        r"""
        The ``n``-th coboundary group of the algebra.

        This is a vector space over the base field `F`, and it is
        returned as a subspace of the vector space `F^d`, where the
        ``n``-th homogeneous component has dimension `d`.

        INPUT:

        - ``n`` -- degree

        EXAMPLES::

            sage: A.<x,y,z> = GradedCommutativeAlgebra(QQ, degrees=(1,1,2))
            sage: d = A.differential({z: x*z})
            sage: d.coboundaries(2)
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: d.coboundaries(3)
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
        """
        A = self.domain()
        F = A.base_ring()
        if n == 0:
            return VectorSpace(F, 0)
        if n == 1:
            return VectorSpace(F, 0)
        M = self.differential_matrix(n-1)
        V0 = VectorSpace(F, M.nrows())
        V1 = VectorSpace(F, M.ncols())
        mor = V0.Hom(V1)(M)
        return mor.image()

    def cocycles(self, n):
        r"""
        The ``n``-th cocycle group of the algebra.

        This is a vector space over the base field `F`, and it is
        returned as a subspace of the vector space `F^d`, where the
        ``n``-th homogeneous component has dimension `d`.

        INPUT:

        - ``n`` -- degree

        EXAMPLES::

            sage: A.<x,y,z> = GradedCommutativeAlgebra(QQ, degrees=(1,1,2))
            sage: d = A.differential({z: x*z})
            sage: d.cocycles(2)
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
        """
        A = self.domain()
        F = A.base_ring()
        if n == 0:
            return VectorSpace(F, 1)
        M = self.differential_matrix(n)
        V0 = VectorSpace(F, M.nrows())
        V1 = VectorSpace(F, M.ncols())
        mor = V0.Hom(V1)(M)
        return mor.kernel()

    def cohomology_raw(self, n):
        r"""
        The ``n``-th cohomology group of ``self``.

        This is a vector space over the base ring, and it is returned
        as the quotient cocycles/coboundaries.

        INPUT:

        - ``n`` -- degree

        .. SEEALSO::

            :meth:`cohomology`

        EXAMPLES::

            sage: A.<x,y,z,t> = GradedCommutativeAlgebra(QQ, degrees=(2,3,2,4))
            sage: d = A.differential({t: x*y, x: y, z: y})
            sage: d.cohomology_raw(4)
            Vector space quotient V/W of dimension 2 over Rational Field where
            V: Vector space of degree 4 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0    0 -1/2]
            [   0    1   -2    1]
            W: Vector space of degree 4 and dimension 0 over Rational Field
            Basis matrix:
            []

        Compare to :meth:`cohomology`::

            sage: d.cohomology(4)
            Free module generated by {[-1/2*x^2 + t], [x^2 - 2*x*z + z^2]} over Rational Field
        """
        return self.cocycles(n).quotient(self.coboundaries(n))

    def cohomology(self, n):
        r"""
        The ``n``-th cohomology group of ``self``.

        This is a vector space over the base ring, defined as the
        quotient cocycles/coboundaries. The elements of the quotient
        are lifted to the vector space of cocycles, and this is
        described in terms of those lifts.

        INPUT:

        - ``n`` -- degree

        .. SEEALSO::

            :meth:`cohomology_raw`

        EXAMPLES::

            sage: A.<a,b,c,d,e> = GradedCommutativeAlgebra(QQ, degrees=(1,1,1,1,1))
            sage: d = A.differential({d: a*b, e: b*c})
            sage: d.cohomology(2)
            Free module generated by {[c*e], [c*d - a*e], [b*e], [b*d], [a*d], [a*c]} over Rational Field

        Compare to :meth:`cohomology_raw`::

            sage: d.cohomology_raw(2)
            Vector space quotient V/W of dimension 6 over Rational Field where
            V: Vector space of degree 10 and dimension 8 over Rational Field
            Basis matrix:
            [ 0  1  0  0  0  0  0  0  0  0]
            [ 0  0  1  0  0  0 -1  0  0  0]
            [ 0  0  0  1  0  0  0  0  0  0]
            [ 0  0  0  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  1  0  0  0  0]
            [ 0  0  0  0  0  0  0  1  0  0]
            [ 0  0  0  0  0  0  0  0  1  0]
            [ 0  0  0  0  0  0  0  0  0  1]
            W: Vector space of degree 10 and dimension 2 over Rational Field
            Basis matrix:
            [0 0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 1]
        """
        H = self.cohomology_raw(n)
        H_basis_raw = [H.lift(H.basis()[i]) for i in range(H.dimension())]
        A = self.domain()
        B = A.basis(n)
        H_basis = [sum([c*b for (c,b) in zip(coeffs, B)]) for coeffs in H_basis_raw]
        # Put brackets around classes.
        H_basis_brackets = [CohomologyClass(b) for b in H_basis]
        return CombinatorialFreeModule(A.base_ring(), H_basis_brackets)
