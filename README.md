nupyck
======
`nupyck` is python frontend to [NUPACK](http://nupack.org), the
nucleic acid package. It has been designed to be familiar
to users of NUPACK, while also presenting a more powerful interactive interface.

Instead of calling NUPACK's executables via the shell, `nupyck` includes
a version of NUPACK compiled as a shared library, which is called directly
from python using `ctypes`.

### Usage
`nupyck` currently exposes four "applications" from the NUPACK suite:
`pfunc`, `pairs`, `concentrations`, and `mfe`. Each of these takes input
similar to their NUPACK counterparts.

#### Options
Common to all four applications is an `Options` object. The defaults are
the same as NUPACK's:

```python
import nupyck

opts = nupyck.Options(
    material=nupyck.RNA,        # | nupyck.DNA
    na=1.0,                     # sodium concentration
    mg=0.0,                     # magnesium concentration
    pseudo=False,               # calculate pseudoknots (only for RNA)
    dangles=nupyck.SOME_DANGLES # | nupyck.NO_DANGLES | nupyck.ALL_DANGLES
)

# or equivalently, because all of the above are defaults:
opts = nupyck.Options()

# or, just set one option:
opts = nupyck.Options(material = nupyck.DNA)
```

#### `pfunc` application
The input's to nupyck's `pfunc` application are the same as to NUPACK's:
a list of sequences, a 1-indexed permutation of those sequences, the
temperature, and other options. For example, consider the following
example from the NUPACK documentation:

```python
import nupyck

sequences = [
    "AGTCTAGGATTCGGCGTGGGTTAA",
    "TTAACCCACGCCGAATCCTAGACTCAAAGTAGTCTAGGATTCGGCGTG",
    "AGTCTAGGATTCGGCGTGGGTTAACACGCCGAATCCTAGACTACTTTG"
]

permutation = [1, 2, 2, 3]

options = nupyck.Options(material = nupyck.DNA)

result = nupyck.pfunc(sequences, permutation, temp = 23, options = options)

print result
```

The result is a python dictionary containing entries for `energy` and `pfunc`,
so the above code produces the following output:
```python
{'energy': -121.97280668378544, 'pfunc': 1.02444390480969e+90}
```
