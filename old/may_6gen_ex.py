
from Nakamura import *
from sage.all import *

def main():
    R = PolynomialRing(GF(2), 'h10,h11,h12,h20,h21,h30')
    gens = [[[(1,0)]],[[(1,1)]],[[(1,2)]],[[(2,0)]],[[(2,1)]],[[(3,0)]]]
    D = Differential(R, gens)
    D.turn_page()
    D.turn_page()

if __name__ == "__main__":
    main()
