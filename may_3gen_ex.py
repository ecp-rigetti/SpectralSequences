from Nakamura import *
from sage.all import *

def main():
    R = PolynomialRing(GF(2), 'h10,h11,h20')
    gens = [[[(1,0)]],[[(1,1)]],[[(2,0)]]]
    D = Differential(R, gens)
    D.turn_page()
    D.turn_page()

if __name__ == "__main__":
    main()

"""
For some reason the last relation (h11^3) keeps adding to relations

"""
