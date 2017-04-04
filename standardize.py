import pandas
import glob
import numpy as np

def read_counts(fn):
    data = pandas.read_csv(fn, sep= '\t', header = 0)
    names  = data[["ID"]].as_matrix().reshape(-1)
    seq    = data[["Sequence"]].as_matrix().reshape(-1)
    counts = data[["Counts"]].as_matrix().reshape(-1).astype('str')
    return names, seq, counts


# Read in pilot design names...
pilot_names = {}
with open('raw/coords_PilotDesign_hg19.txt') as fp:
    for line in fp:
        name, chrom, start, end = line.strip().split()
        pilot_names[name] = "{}:{}-{}:+:{}".format(chrom, start, end, name)




# GSM1831757_HepG2_ScaleUpDesign1_minP_mRNA_Rep1.counts.txt.gz
cells = ["HepG2", "K562"]
designs = ["PilotDesign", "ScaleUpDesign1", "ScaleUpDesign2"]
promoters = ["minP", "SV40P"]

for cell in cells:
    for promoter in promoters:
        for design in designs:
            # Pilot only done with SV40P
            if design == designs[0] and promoter == promoters[0]: continue

            # Read in counts.
            rna1 = glob.glob("raw/*{}_{}_{}_mRNA_Rep1*".format(cell, design, promoter))
            rna2 = glob.glob("raw/*{}_{}_{}_mRNA_Rep2*".format(cell, design, promoter))
            dna = glob.glob("raw/*{}_{}_Plasmid*".format(design, promoter))
            print rna1, rna2, dna
            assert len(rna1) == len(rna2) and 0 < len(dna) < 3
            rna1 = rna1[0]
            rna2 = rna2[0]
            dna1 = dna[0]
            dna2 = dna[0] if len(dna) == 1 else dna[1]
            n1, s1, c1 = read_counts(rna1)
            n2, s2, c2 = read_counts(rna2)
            nd1, sd1, cd1 = read_counts(dna1)
            nd2, sd2, cd2 = read_counts(dna2)
            assert np.all(s1 == s2) and np.all(s2 == sd1) and np.all(sd1 == sd2)
            assert np.all(n1 == n2) and np.all(n2 == nd1) and np.all(nd1 == nd2)
            assert c1.shape == c2.shape == cd1.shape == cd2.shape
            seqs, names = s1, n1

            # Infer genomic position from name.
            std_names = []
            for name in names:
                if design == designs[0]:
                    std_names += [pilot_names[name]]
                else:
                    ofs, chrom, base = name.split('_')[3:]
                    ofs, base = int(ofs), int(base)
                    assert chrom[:3] == 'chr'
                    idx = base + ofs * 5
                    std_names += ["{}:{}-{}:+:{}".format(chrom, idx, idx + 145, name)]
                
            # Write counts file.
            with open("standard/{}_{}_{}.counts".format(cell, design, promoter), 'w') as fp:
                fp.write('\t'.join(["NAME", "DNA", "RNA"]) + '\n')
                for name, d1, d2, r1, r2 in zip(std_names, cd1, cd2, c1, c2):
                    d = d1 + ':' +  d2
                    r = r1 + ':' +  r2
                    fp.write('\t'.join([name, d, r]) + '\n')
            
            # Write fasta file.
            with open("standard/{}_{}_{}.fa".format(cell, design, promoter), 'w') as fp:
                for name, seq in zip(std_names, seqs):
                    fp.write('>' + name + '\n')
                    fp.write(seq + '\n')
