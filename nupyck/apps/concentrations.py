from .. import core
import pfunc
import numpy as np
import itertools

def melting_temp(template, primer, t_conc, p_conc, t_lo, t_hi, na=1.0, mg=0.0,
        eps = 0.02):
    def get_frac(temp):
        return tp_fraction(template, primer, t_conc, p_conc, temp, na, mg)

    # base case: not gonna get better
    frac = get_frac(t_lo)
    if frac <= 0.5:
        return t_lo

    # base case: not gonna get worse
    frac = get_frac(t_hi)
    if frac >= 0.5:
        return t_hi

    converged = False
    while not converged:
        temp = t_lo + (t_hi - t_lo) / 2.
        frac = get_frac(temp)

        # search lower
        if frac < 0.5 - eps:
            t_hi = temp

        # search higher
        elif frac > 0.5 + eps:
            t_lo = temp

        else:
            converged = True

    return temp

def _convert_perms_to_A(perms):
    numTotal = len(perms)
    numSS = max(map(max, perms))

    A = np.zeros((numSS, numTotal), dtype=int)

    for p, perm in enumerate(perms):
        for ss in perm:
            A[ss-1][p] += 1

    return A

def tp_fraction(template, primer, t_conc, p_conc, temp, na=1.0, mg=0.0):
    """Return the fraction of template converted to template/primer duplex"""

    x0 = np.array([t_conc, p_conc])

    G = np.array(
            [ pfunc.pfunc(
                [ template, primer ],
                temp = temp, na = na, mg = mg,
                material = core.DNA,
                perm = p
              )['energy'] for p in [1], [2], [1,1], [1,2], [2,2]
            ]
          )

    A  = np.array([[1, 0, 2, 1, 0], [0, 1, 0, 1, 2]])

    x = calc_conc(x0, G, A, temp)
    tp_conc = x[3]
    frac = tp_conc / t_conc

    return frac

def concentrations(
        species,
        x0,
        max_complex_size,
        temp=37,
        names=None,
        options=core.Options()):

    perms = itertools.chain(
        *[itertools.combinations_with_replacement(range(1,len(species)+1), n)
            for n in range(1, max_complex_size+1)
            ]
        )
    perms = list(perms)

    G = np.array(
            [ pfunc.pfunc(species, permutation=p, temp=temp, options=options)['energy']
                for p in perms
                ]
            )

    A = _convert_perms_to_A(perms)

    x = calc_conc(np.array(x0), G, A, temp)

    if names is None:
        names = perms
    else:
        names = ["-".join([names[p-1] for p in perm]) for perm in perms]

    return {
        "concentrations": dict(zip(names, x)),
        "energies": dict(zip(names, G))
    }
    return dict(zip(names, x))


def calc_conc(x0, G, A, temp):

    numSS, numTotal = A.shape
    kT = (273.15 + temp) * core.kB
    waterDensity = core.nupack.WaterDensity(core.c_double(temp))

    x = core.c_double_array([-1] * numTotal)
    A            = (core.POINTER(core.c_int) * numSS)(
                       *[ core.c_array(a) for a in A ]
                   )
    G            = core.c_double_array(G / kT)
    x0           = core.c_double_array(x0 / waterDensity)
    numSS        = core.c_int(numSS)
    numTotal     = core.c_int(numTotal)
    maxIters     = core.c_int(10000)
    tol          = core.c_double(1e-7)
    deltaBar     = core.c_double(1000)
    eta          = core.c_double(0.125)
    kT           = core.c_double(kT)
    maxNoStep    = core.c_int(50)
    maxTrial     = core.c_int(100000)
    perturbScale = core.c_double(100)
    quiet        = core.c_int(1)
    writeLogFile = core.c_int(0)
    logFile      = core.c_char_p(None)
    h20Density   = core.c_double(waterDensity)
    seed         = core.c_int(0)

    converged = core.nupack.CalcConc(
      x, A, G, x0, numSS, numTotal, maxIters, tol, deltaBar,
      eta, kT, maxNoStep, maxTrial, perturbScale, quiet, writeLogFile,
      logFile, h20Density, seed
    )

    if not converged:
      raise RuntimeError("concentration calculation did not converge")

    x = np.array(x[:]) * waterDensity
    return x
