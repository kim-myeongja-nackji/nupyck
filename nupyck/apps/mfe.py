from .. import core


class oneDnaStruct(core.Structure):
    _fields_ = ([
      ("theStruct", core.POINTER(core.c_int)),
      ("error", core.c_longdouble),
      ("correctedEnergy", core.c_longdouble),
      ("slength", core.c_int)
    ])


class dnaStructures(core.Structure):
    _fields_ = ([
      ("validStructs", core.POINTER(oneDnaStruct)),
      ("nStructs", core.c_int),
      ("nAlloc", core.c_int),
      ("seqlength", core.c_int),
      ("minError", core.c_longdouble)
    ])

def dotParen(structure, sequence):

    dP = ['.'] * structure.slength
    for i in range(len(dP)):
        if structure.theStruct[i] > -1 and dP[i] == '.':
            dP[i] = '('
            dP[structure.theStruct[i]] = ')'

    result = list(sequence)
    dP_index = 0
    for i in range(len(sequence)):
        if sequence[i] != '+':
            result[i] = dP[dP_index]
            dP_index += 1

    return "".join(result)


def mfe(sequences,
        permutation,
        temp=37,
        degenerate=False,
        options=core.Options()):

    sequence, symmetry = options.joinSequences(sequences, permutation)

    seq_as_ints = core.seqToInts(sequence)
    seq_len = len(sequence) - sequence.count('+')

    mfe_structs = dnaStructures(
        core.POINTER(oneDnaStruct)(),
        core.c_int(0),
        core.c_int(0),
        core.c_int(0),
        core.c_longdouble(100000) # "infinity"
    )

    complexity = 5 if options.pseudo else 3

    only_one = 0 if degenerate else 1

    core.nupack.mfeFullWithSym(
        core.c_array(seq_as_ints),     # int inputSeq[]
        core.c_int(seq_len),           # int seqLen
        core.byref(mfe_structs),       # dnaStructures *mfeStructures
        core.c_int(complexity),        # int complexity
        core.c_int(options.material),  # int naType
        core.c_int(options.dangles),   # int dangles
        core.c_longdouble(temp),       # double temperature
        core.c_int(symmetry),          # int symmetry
        core.c_int(only_one),          # int onlyOne
        core.c_longdouble(options.na), # double sodiumconc
        core.c_longdouble(options.mg), # double magnesiumconc
        core.c_longdouble(0)           # int uselongsalt
    )

    results = [
        { "structure": dotParen(structure, sequence),
          "energy"   : structure.correctedEnergy
        }
        for structure in mfe_structs.validStructs[:mfe_structs.nStructs]
    ]

    return results
