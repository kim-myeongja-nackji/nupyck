from .. import core
import pfunc
import numpy as np

def melting_temp(template, primer, t_conc, p_conc, trange, na=1.0, mg=0.0):
    maxfrac = 0
    for tidx,temp in enumerate(trange):
        frac = tp_fraction(template, primer, t_conc, p_conc, temp, na, mg)
        if frac < 0.5:
            return trange[max(tidx-1, 0)]
    return trange[-1]

def tp_fraction(template, primer, t_conc, p_conc, temp, na=1.0, mg=0.0):
  """Return the fraction of template converted to template/primer duplex"""

  t  = template
  p  = primer
  tt = "+".join((template, template))
  tp = "+".join((template, primer))
  pp = "+".join((primer, primer))

  kT = (temp + 273.15) * core.kB

  G = [ pfunc.single(
          seq, temp, material = core.DNA, na = na, mg = mg, symmetry = sym
        )['energy'] / kT
        for seq, sym in zip([ t, p, tt, tp, pp ], [ 1, 1, 2, 1, 2 ])
      ]

  waterDensity = core.nupack.WaterDensity(core.c_double(temp))

  # result array
  x = core.c_double_array([-1,-1,-1,-1,-1])

  A = (core.POINTER(core.c_int) * 2)(
    core.c_array([1, 0, 2, 1, 0]),
    core.c_array([0, 1, 0, 1, 2])
  )
  G            = core.c_double_array(G)
  x0           = core.c_double_array([ t_conc/waterDensity, p_conc/waterDensity ])
  numSS        = core.c_int(2)
  numTotal     = core.c_int(5)
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

  # convert back to molarity
  x = np.array(x) * waterDensity

  # final [template-primer duplex] vs initial [template]
  conc_ratio = x[3] / t_conc
  return conc_ratio

