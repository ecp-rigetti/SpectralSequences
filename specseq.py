from sage.all import *
from cdga import *

class Spectral_Sequence_char2:

    def __init__(self, gens, degrees, diff):
        self.page = 1
        self.CDGA = CDGA_char2(gens, diff, degrees)
        self.original_ring = self.CDGA.R


    def turn_page(self, new_diffs):
        cycles, bds = self.CDGA.reduction_step()

        self.page += 1
