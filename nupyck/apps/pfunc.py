from .. import core
import math
import numpy as np

calcVPi = core.nupack.calculateVPi
calcVPi.restype = core.c_int

def pfunc(sequence, temp=37,
        material=core.RNA,
        na=1.0, mg=0.0,
        pseudo=False,
        dangles=core.SOME_DANGLES,
        calc_pairs = False,
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

    if calc_pairs:
        if pseudo:
            raise NotImplementedError

        mat_dim = len(sequence) - sequence.count("+")
        mat_size = mat_dim * (mat_dim + 1)
        pair_pr = (core.c_longdouble * mat_size)()
        core.pair_pr.contents = pair_pr

    seq_as_ints = core.seqToInts(sequence)
    complexity = 5 if pseudo else 3

    pf = core.nupack.pfuncFullWithSym(
      core.c_array(seq_as_ints), # inputSeq
      core.c_int(complexity),    # complexity
      core.c_int(material),      # naType
      core.c_int(dangles),       # dangles
      core.c_longdouble(temp),   # temperature
      core.c_int(calc_pairs),    # calcPairs
      core.c_int(symmetry),      # permSymmetry
      core.c_longdouble(na),     # sodiumconc
      core.c_longdouble(mg),     # magnesiumconc
      core.c_int(0)              # uselongsalt
    )

    energy = -core.kB * (273.15 + temp) * math.log(max(pf,1))

    if calc_pairs:
        ppairs = np.array(pair_pr[:mat_size]).reshape(mat_dim, mat_dim + 1)
        return {'energy' : energy, 'pfunc' : pf, 'ppairs' : ppairs}

    else:
        return {'energy' : energy, 'pfunc' : pf}

def pairs(seq, **kwargs):
    kwargs['calc_pairs'] = True
    return pfunc(seq, **kwargs)
