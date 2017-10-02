from .. import core
import math

calcVPi = core.nupack.calculateVPi
calcVPi.restype = core.c_int

def pfunc(sequence, temp=37,
        material=core.RNA,
        na=1.0, mg=0.0,
        pseudo=False,
        dangles=core.SOME_DANGLES,
        perm = None):

    if perm is None:
        symmetry = 1

    else:
        if type(sequence) is not list:
            raise ValueError("permutation given without list of sequences")

        sequence = "+".join(
            [ sequence[p-1] for p in perm ]
        )

        symmetry = calcVPi(core.c_array(perm), core.c_int(len(perm)))

    if pseudo and material == core.DNA:
        raise ValueError("pseudoknot option valid for RNA only")

    if pseudo and '+' in sequence:
        raise ValueError("pseudoknot option valid only for single strands")

    if material == core.RNA and (na != 1.0 or mg != 0.0):
        raise ValueError("salt corrections unavailable for RNA")

    seq_as_ints = core.seqToInts(sequence)
    complexity = 5 if pseudo else 3

    pf = core.nupack.pfuncFullWithSym(
      core.c_array(seq_as_ints), # inputSeq
      core.c_int(complexity),    # complexity
      core.c_int(material),      # naType
      core.c_int(dangles),       # dangles
      core.c_longdouble(temp),   # temperature
      core.c_int(False),         # calcPairs
      core.c_int(symmetry),      # permSymmetry
      core.c_longdouble(na),     # sodiumconc
      core.c_longdouble(mg),     # magnesiumconc
      core.c_int(0)              # uselongsalt
    )

    energy = -core.kB * (273.15 + temp) * math.log(max(pf,1))
    return {'energy' : energy, 'pfunc' : pf}

