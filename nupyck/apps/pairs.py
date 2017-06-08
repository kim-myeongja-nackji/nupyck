from .. import core
import pfunc
import numpy as np

def single(seq, temp, *args, **kwargs):

    if kwargs.get('pseudo', False):
        raise NotImplementedError(
            "pair probability with pseudoknots not implemented yet."
        )

    seq_len = len(seq)
    mat_size = seq_len * (seq_len + 1)
    pair_pr = core.POINTER(core.c_longdouble).in_dll(core.nupack, "pairPr")
    pair_pr.contents = (core.c_longdouble * mat_size)()

    pfunc.single(seq, temp, calc_pairs = True, *args, **kwargs)

    return np.array(pair_pr[:mat_size]).reshape(seq_len, seq_len + 1)

def multi(seqs, temp, *args, **kwargs):

    # TODO
    raise NotImplementedError("TODO")
