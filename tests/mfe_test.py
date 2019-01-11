import context
import nupyck

sequences = [
    "AGTCTAGGATTCGGCGTGGGTTAA",
    "TTAACCCACGCCGAATCCTAGACTCAAAGTAGTCTAGGATTCGGCGTG",
    "AGTCTAGGATTCGGCGTGGGTTAACACGCCGAATCCTAGACTACTTTG"
]

options = nupyck.Options(material = nupyck.DNA)

res = nupyck.mfe(sequences, [1, 2, 2, 3], 23, options = options)

print res
