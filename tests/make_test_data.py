import numpy as np
import sys, os, subprocess, tempfile, json, argparse

nupack_dir = "../lib/nupack"
os.environ['NUPACKHOME'] = nupack_dir
os.environ['PATH'] = os.environ['PATH'] + ':%s' % (nupack_dir + '/bin')

def randseq(n):
    return "".join(np.random.choice(list("ATCG"), n))

def pfunc(num_samples):

    jobs = []

    for _ in range(num_samples):

        n_seqs = np.random.randint(1, 4)
        seq_lengths = np.random.randint(20, 40, size=n_seqs)
        sequences = [ randseq(n) for n in seq_lengths ]

        perm_size = np.random.randint(1, n_seqs + 3)
        perm = list(np.random.randint(1, n_seqs+1, size=perm_size))

        temp = np.random.uniform(0, 100)

        material = np.random.choice(["rna", "dna"])

        job_input = (
            "%d\n"
            "%s\n"
            "%s\n"
        ) % (
            n_seqs,
            "\n".join(sequences),
            " ".join(map(str, perm))
        )

        with tempfile.NamedTemporaryFile(suffix='.in') as input_file:

            input_file.write(job_input)
            input_file.flush()

            job_command = "pfunc -material %s -T %f -multi %s" % (
                material, temp, input_file.name.replace('.in', '')
            )

            out = subprocess.check_output(job_command.split(' '))

        energy, pf = [
            float(line)
            for line in out.strip().split('\n')
            if not line.startswith('%')
        ]

        jobs.append({
            "sequences": sequences,
            "perm": perm,
            "temp": temp,
            "material": material,
            "energy": energy,
            "pfunc": pf
        })

    return jobs



if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("app", choices=["pfunc"])
    parser.add_argument("--num-samples", type=int, default=100, metavar="N")
    parser.add_argument("--seed", type=int)

    args = parser.parse_args()

    if args.seed:
        np.random.seed(args.seed)

    jobs = eval(args.app)(args.num_samples)

    json.dump(jobs, sys.stdout, indent=1)

