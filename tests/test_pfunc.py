import context
import nupyck
import json
import sys

jobs = json.load(sys.stdin)

print "Checking %d inputs" % len(jobs)
all_ok = True

for job in jobs:
    material = nupyck.DNA if job['material'] == 'dna' else nupyck.RNA
    options = nupyck.Options(material = material)
    res = nupyck.pfunc(job['sequences'], job['perm'], job['temp'], options = options)
    pfunc_err = abs(1.0 - res['pfunc'] / job['pfunc'])
    energy_err = abs(res['energy'] - job['energy'])

    if pfunc_err > 1e-5 or energy_err > 1e-3:
        all_ok = False
        print "Check failed for this input:"
        print job
        print "Gave result:"
        print res

if all_ok:
    print "All tests passed."
