from .. import core
import pfunc
import numpy as np

def _call_pf(seq, temp, mat_dim, *args, **kwargs):
    mat_size = mat_dim * (mat_dim + 1)

    pair_pr = core.POINTER(core.c_longdouble).in_dll(core.nupack, "pairPr")
    pair_pr.contents = (core.c_longdouble * mat_size)()

    pf = pfunc.single(seq, temp, calc_pairs = True, *args, **kwargs)

    ppairs = np.array(pair_pr[:mat_size]).reshape(mat_dim, mat_dim + 1)

    return { 'ppairs' : ppairs, 'pfunc' : pf }


def single(seq, temp, *args, **kwargs):

    if kwargs.get('pseudo', False):
        raise NotImplementedError(
            "pair probability with pseudoknots not implemented yet."
        )

    mat_dim = len(seq)
    return _call_pf(seq, temp, mat_dim, *args, **kwargs)


def multi(seqs, perm, temp, *args, **kwargs):

    if kwargs.get('pseudo', False):
        raise ValueError("pseudoknot option valid only for single strands.")

    seq_perm = [ seqs[p-1] for p in perm ]
    mat_dim = sum(map(len, seq_perm))
    seq = "+".join(seq_perm)

    return _call_pf(seq, temp, mat_dim, *args, **kwargs)
