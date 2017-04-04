"""
Accepts an experiment basename given as argv[1].
Calculates normalized activity scores for each
entry as log(RNA+1) - log(DNA+1).
Writes data to argv[2] with format
NAME\tRep1_m1,Rep1_m2,...:Rep2_m1,Rep2_m2,...

TODO: either expose the DNA count or
estimate the variance of each measurement here.
"""
from math import log
THRESH = 20

def normalize(counts, out):
    with open(counts) as fp, open(out, 'w') as out_fp:
        fp.readline()
        out_fp.write('NAME\tACTIVITY\n')
        for line in fp:
            name, dna, rna = line.strip().split()
            norm = []
            for exp_dna, exp_rna in zip(dna.split(':'), rna.split(':')):
                exp_norm = []
                for rep_dna, rep_rna in zip(exp_dna.split(','), exp_rna.split(',')):
                    if int(rep_dna) > 20:
                        exp_norm += [log(float(rep_rna) + 1.0) / log(float(rep_dna) + 1.0)]
                    else:
                        exp_norm += [float('inf')]
                norm += [exp_norm]
            out_fp.write(name + '\t' + ':'.join(map(lambda x: ','.join(map(str, x)), norm)) + '\n')

if __name__ == '__main__':
    import sys
    normalize(*sys.argv[1:])
