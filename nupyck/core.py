import numpy as np
import os, math

from ctypes import *

package_dir = os.path.dirname(__file__)
nupack_dir  = os.path.join(package_dir, '../lib/nupack')
os.environ['NUPACKHOME'] = nupack_dir

nupack = cdll.LoadLibrary(os.path.join(nupack_dir, 'nupack.so'))
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

class oneDnaStruct(Structure):
  _fields_ = ([
    ("theStruct", POINTER(c_int)),
    ("error", c_longdouble),
    ("correctedEnergy", c_longdouble),
    ("slength", c_int)
  ])

class dnaStructures(Structure):
  _fields_ = ([
    ("validStructs", POINTER(oneDnaStruct)),
    ("nStructs", c_int),
    ("nAlloc", c_int),
    ("seqlength", c_int),
    ("minError", c_longdouble)
  ])



pairs = dict(zip("GATC", "CTAG"))
def pairof(b): return pairs[b]

def tp_end_pairfrac(template, primer, t_conc, p_conc, temp, na=1.0, mg=0.0):
  lp = len(primer)

  if template[-lp] != pairof(primer[-1]):
    return 0.0

  lt = len(template)
  ppair_idx1 = lt - lp
  ppair_idx2 = lt + lp - 1

  t  = seqToInts(template)
  p  = seqToInts(primer)
  tt = seqToInts("+".join((template, template)))
  tp = seqToInts("+".join((template, primer)))
  pp = seqToInts("+".join((primer, primer)))

  ppair, tp_pf = ppair_single(
    tp, temp,
    ppair_idx1, ppair_idx2,
    na = na, mg = mg,
    get_pf = True
  )
  tp_energy = -math.log(max(tp_pf, 1))

  get_energy = lambda s: -math.log(max(pfunc(s, temp, na = na, mg = mg), 1))
  G = ([
    get_energy(t),
    get_energy(p),
    get_energy(tt),
    tp_energy,
    get_energy(pp)
  ])

  waterDensity = nupack.WaterDensity(c_double(temp))

  # result array
  x = c_double_array([-1,-1,-1,-1,-1])

  A = (POINTER(c_int) * 2)(
    c_array([1, 0, 2, 1, 0]),
    c_array([0, 1, 0, 1, 2])
  )
  G            = c_double_array(G)
  x0           = c_double_array([ t_conc/waterDensity, p_conc/waterDensity ])
  numSS        = c_int(2)
  numTotal     = c_int(5)
  maxIters     = c_int(10000)
  tol          = c_double(1e-7)
  deltaBar     = c_double(1000)
  eta          = c_double(0.125)
  kT           = c_double((temp + 273.15) * kB)
  maxNoStep    = c_int(50)
  maxTrial     = c_int(100000)
  perturbScale = c_double(100)
  quiet        = c_int(1)
  writeLogFile = c_int(0)
  logFile      = c_char_p(None)
  h20Density   = c_double(waterDensity)
  seed         = c_int(0)

  converged = nupack.CalcConc(
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

  pfrac = ppair * conc_ratio

  return pfrac


def mfe(seq_as_ints, temp, na = 1.0, mg = 0.0):

  mfeStructs = dnaStructures(
    POINTER(oneDnaStruct)(), 
    c_int(0),
    c_int(0),
    c_int(0),
    c_longdouble(100000)
  )

  nupack.mfeFullWithSym(
    c_array(seq_as_ints),
    c_int(len(seq_as_ints)-1),
    byref(mfeStructs),
    c_int(3),
    c_int(0),
    c_int(1),
    c_longdouble(temp),
    c_int(1),
    c_int(1),
    c_longdouble(na),
    c_longdouble(mg),
    c_int(0)
  )

  return mfeStructs.validStructs[:mfeStructs.nStructs]

def mfe_prob(seq_as_ints, temp):
  pf = pfunc(seq_as_ints, temp)
  energy = mfe(seq_as_ints, temp)[0].correctedEnergy

  return math.exp(-energy / (kB * (temp + 273.15))) / pf

def display_struct(seq_as_ints, dnaStruct):
  nickPosition = seq_as_ints.index(15)

  slength   = dnaStruct.slength
  theStruct = dnaStruct.theStruct

  result = ['.'] * slength
  for i in range(slength):
    if theStruct[i] > -1 and result[i] == '.':
      result[i] = '('
      result[theStruct[i]] = ')'

  result.insert(nickPosition, '+')

  print intsToSeq(seq_as_ints)
  print "".join(result)

