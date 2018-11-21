from .. import core
import math
import numpy as np

def pfunc(sequences,
          permutation,
          temp=37,
          calc_pairs=False,
          options=core.Options()):

    sequence, symmetry = options.joinSequences(sequences, permutation)

    if calc_pairs:
        if options.pseudo:
            raise NotImplementedError(
                "pair calculation for pseudoknots is not currently "+
                "implemented in nupyck"
            )

        mat_dim = len(sequence) - sequence.count("+")
        mat_size = mat_dim * (mat_dim + 1)
        pair_pr = (core.c_longdouble * mat_size)()
        core.pair_pr.contents = pair_pr

    seq_as_ints = core.seqToInts(sequence)
    complexity = 5 if options.pseudo else 3

    pf = core.nupack.pfuncFullWithSym(
      core.c_array(seq_as_ints),     # inputSeq
      core.c_int(complexity),        # complexity
      core.c_int(options.material),  # naType
      core.c_int(options.dangles),   # dangles
      core.c_longdouble(temp),       # temperature
      core.c_int(calc_pairs),        # calcPairs
      core.c_int(symmetry),          # permSymmetry
      core.c_longdouble(options.na), # sodiumconc
      core.c_longdouble(options.mg), # magnesiumconc
      core.c_int(0)                  # uselongsalt
    )

    energy = -core.kB * (273.15 + temp) * math.log(max(pf,1))

    if calc_pairs:
        ppairs = np.array(pair_pr[:mat_size]).reshape(mat_dim, mat_dim + 1)
        return {'energy' : energy, 'pfunc' : pf, 'ppairs' : ppairs}

    else:
        return {'energy' : energy, 'pfunc' : pf}

def pairs(*args, **kwargs):
    kwargs['calc_pairs'] = True
    return pfunc(*args, **kwargs)
