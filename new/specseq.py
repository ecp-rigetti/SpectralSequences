from sage.all import *
from cdga import *
from MayDifferential import *

class Spectral_Sequence_char2:

    """
    initialize a spectral sequence class
    diff should be a differential class that can compute differentials
    at different pages
    """
    def __init__(self, gens, degrees, diff, relts):
        self.page = 1
        self.n = 0
        #self.Differential = diff
        #first_diff = self.Differential.get_diff()
        self.relts = relts
        self.CDGA = CDGA_char2(gens, diff, degrees, relts)
        self.original_ring = self.CDGA.R
        self.dict = {}

    def new_varname(self):
        self.n += 1
        return 'var' + str(self.n)

    def turn_page(self):
        cycles, bds = self.CDGA.reduction_step()
        print cycles, bds
        degrees = map(lambda c: c.degree(), cycles)
        self.page += 1
        new_diff = []
        for c in cycles:
            k = input('input the {}th differential for {}\n'.format(self.page, c))
            new_diff.append(k)
        new_gens = []
        for c in cycles:
            if c.is_generator():
                new_gens.append(str(c))
            else:
                new_var = self.new_varname()
                new_gens.append(new_var)
                self.dict[new_var] = str(c)
        diff = {}
        for i in range(len(new_gens)):
            diff[new_gens[i]] = new_diff[i]
        new_relts = self.relts + map(lambda b: str(b), bds)
        self.CDGA = CDGA_char2(new_gens, diff, degrees, new_relts)


def main():
    A = Spectral_Sequence_char2('h0,h1,h2,h20,h21,h30', degrees=None, diff={'h20': 'h0*h1', 'h21': 'h1*h2', 'h30': 'h21*h0 + h2*h20'}, relts=[])
    A.turn_page()
    # B = Spectral_Sequence_char2('h0,h1,b20', degrees=(1,1,2), diff={'b20':'h1**3'},relts=['h0*h1'])
    # B.turn_page()



if __name__ == "__main__":
    main()
