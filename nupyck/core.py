import numpy as np
import os, math

from ctypes import *

package_dir = os.path.dirname(__file__)
os.environ['NUPACKHOME'] = package_dir

nupack = cdll.LoadLibrary(os.path.join(package_dir, 'nupack.so'))
nupack.pfuncFull.restype = c_longdouble
nupack.pfuncFullWithSym.restype = c_longdouble
nupack.pfunc.restype     = c_longdouble
nupack.WaterDensity.restype = c_double

pair_pr = POINTER(c_longdouble).in_dll(nupack, "pairPr")

# boltzmann constant
kB = 0.0019872041

NO_DANGLES = 0
SOME_DANGLES = 1
ALL_DANGLES = 2

DNA = 0
RNA = 1

c_array = lambda l: (c_int * len(l))(*l)
c_double_array = lambda l: (c_double * len(l))(*l)

baseToInt = dict( zip( "ACGTU+" , [1,2,3,4,4,15] ) )
def seqToInts(seq): return [ baseToInt[b] for b in seq ] + [-1]

intToBase = dict( zip( [1,2,3,4,15] , "ACGT+" ) )
def intsToSeq(ints): return "".join([ intToBase[i] for i in ints[:-1] ])

complement = dict( zip("ACGT", "TGCA") )
revcomp = lambda seq: "".join(reversed([complement[b] for b in seq]))


calcVPi = nupack.calculateVPi
calcVPi.restype = c_int


class Options:
    def __init__(
            self,
            material=RNA,
            na=1.0, mg=0.0,
            pseudo=False,
            dangles=SOME_DANGLES):

        self.material = material
        self.na = na
        self.mg = mg
        self.pseudo = pseudo
        self.dangles = dangles

        if pseudo and material == DNA:
            raise ValueError("pseudoknot option valid for RNA only")

        if material == RNA and (na != 1.0 or mg != 0.0):
            raise ValueError("salt corrections unavailable for RNA")

    def joinSequences(self, sequences, permutation):
        if len(permutation) > 1 and self.pseudo:
            raise ValueError("pseudoknot option valid only for single strands")

        sequence = "+".join(
            sequences[p-1] for p in permutation
        )

        symmetry = calcVPi(
            c_array(permutation),
            c_int(len(permutation))
        )

        return sequence, symmetry
