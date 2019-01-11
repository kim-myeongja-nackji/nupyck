import context
import nupyck

sequences = [
    "AGTCTAGGATTCGGCGTGGGTTAA",
    "TTAACCCACGCCGAATCCTAGACTCAAAGTAGTCTAGGATTCGGCGTG",
    "AGTCTAGGATTCGGCGTGGGTTAACACGCCGAATCCTAGACTACTTTG"
]

options = nupyck.Options(material = nupyck.DNA)

res = nupyck.concentrations(
    sequences,
    x0=[1e-6, 1e-6, 1e-6],
    max_complex_size=3
)

print res
