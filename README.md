# Spectral Sequences

This is code from an (ongoing) project begun in Spring 2017 with two goals in mind:

(1) Given a polynomial commutative differential graded algebra over a finite field, compute a minimal basis for its cohomology (i.e., find a minimal generating set for cycles mod boundaries, and a minimal generating set for boundaries).

(2) Using this sequentially, calculate the pages of a spectral sequence in an automated fashion; for this, the toy example of the May spectral sequence was used.

At various points in development, the project has gotten pretty close to successful in both (1) and (2); at the moment, one of the limiting factors is in the final reduction step for the cycles/boundaries generating sets, where we wish to find a minimal multiplicative generating set for the cohomology. The standard division algorithm works in the ideal setting (where a minimal basis corresponds to a set such that a linear combination of its elements can yield any member of the ring), but I'm pretty sure we want a generating set where products and sums of products of generating elements can yield any element of the cohomology.

I'll be continuing to think about this project in the future; in the meantime, if you stumble upon this and (a) it could be of use to you and your research, or (b) you have any insights into how one might enact a reduction to a minimal multiplicative generating set, please feel free to reach out.


## Code:

New code (with a purpose-built "CDGA" class) is found in the new folder

Old code (where everything was kind of jumpled together) is found in the old folder

Both of them work to varying degrees, but they have the same issue with the final reduction/minimal generating set step.


## Dependencies are:

Sage (http://www.sagemath.org/)

(old dependency: Macaulay2 (http://www.math.uiuc.edu/Macaulay2/))

Run python files by invoking python with sage ($ sage -python *.py)

##TODO:

fix up the documentation (old code is pretty well documented, need to update for new code)
