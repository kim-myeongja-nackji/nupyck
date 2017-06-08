from .. import core
import math

def single(sequence, temp,
        material=core.RNA,
        na=1.0, mg=0.0,
        pseudo=False,
        dangles=core.SOME_DANGLES,
        calc_pairs = False):

    if pseudo and material == core.DNA:
        raise ValueError("pseudoknot option valid for RNA only")

    seq_as_ints = core.seqToInts(sequence)
    complexity = 5 if pseudo else 3

    pf = core.nupack.pfuncFull(
      core.c_array(seq_as_ints), # inputSeq
      core.c_int(complexity),    # complexity
      core.c_int(material),      # naType
      core.c_int(dangles),       # dangles
      core.c_longdouble(temp),   # temperature
      core.c_int(calc_pairs),    # calcPairs
      core.c_longdouble(na),     # sodiumconc
      core.c_longdouble(mg),     # magnesiumconc
      core.c_int(0)              # uselongsalt
    )

    energy = -core.kB * (273.15 + temp) * math.log(pf)
    return {'energy' : energy, 'pfunc' : pf}

def multi(sequences, perm, *args, **kwargs):

    if kwargs.get('pseudo', False):
        raise ValueError("pseudoknot option valid only for single strands")

    seq = "+".join(
        [ sequences[p-1] for p in perm ]
    )

    return single(seq, *args, **kwargs)
