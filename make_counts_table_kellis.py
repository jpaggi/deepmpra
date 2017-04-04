import pandas
import glob
import numpy as np

DATA_DIR = '../enhancer_data/kellis/raw'
OUT_DIR = '../enhancer_data/kellis/counts'

cells = ["HepG2", "K562"]
designs = ["ScaleUpDesign1", "ScaleUpDesign2"]
promoters = ["minP", "SV40P"]

def read_counts(fn):
    data = pandas.read_csv(fn, sep= '\t', header = 0)
    names  = data[["ID"]].as_matrix().reshape(-1)
    seq    = data[["Sequence"]].as_matrix().reshape(-1)
    counts = data[["Counts"]].as_matrix().reshape(-1).astype('int')
    return names, seq, counts

def get_counts(cell, promoter, design, names=None, seqs=None):
    # Find relevant files.
    # Note that there is one plasmid library used for each promoter.
    rna1 = glob.glob("{}/*{}_{}_{}_mRNA_Rep1*".format(DATA_DIR, cell, design, promoter))
    rna2 = glob.glob("{}/*{}_{}_{}_mRNA_Rep2*".format(DATA_DIR, cell, design, promoter))
    dna  = glob.glob("{}/*{}_{}_Plasmid*".format(DATA_DIR, design, promoter))
    print rna1, rna2, dna
    assert len(rna1) == len(rna2) == len(dna) == 1

    # Read in names, sequences, and counts.
    rna1_names, rna1_seqs, rna1_counts = read_counts(rna1[0])
    rna2_names, rna2_seqs, rna2_counts = read_counts(rna2[0])
    dna_names,  dna_seqs,  dna_counts  = read_counts(dna[0])

    # ... and assure that they are concordant.
    assert np.all(rna1_seqs == rna2_seqs) and np.all(rna2_seqs == dna_seqs)
    assert np.all(rna1_names == rna2_names) and np.all(rna2_names == dna_names)
    assert rna1_counts.shape == rna2_counts.shape == dna_counts.shape

    # ... and they match the previously read names and sequences.
    if names is not None:
        assert np.all(names == dna_names)
        assert np.all(seqs  == dna_seqs)
    seqs, names = dna_seqs, dna_names
    return names, seqs, rna1_counts, rna2_counts, dna_counts

###########################################################################
######## Read in counts. ##################################################

def read_data():
    names = {}
    seqs = {}
    rna1 = {}
    rna2 = {}
    dna  = {}
    for design in designs:
        names[design] = None
        seqs[design]  = None
        rna1[design] = []
        rna2[design] = []
        dna[design] = []
        header = []
        for promoter in promoters:
            for cell in cells:
                header += ["{}-{}".format(promoter, cell)]
                (names[design], seqs[design],
                 r1, r2, d) = get_counts(cell, promoter, design,
                                         names[design], seqs[design])
                rna1[design] += [r1]
                rna2[design] += [r2]
                dna[design] += [d]
        rna1[design] = np.vstack(rna1[design]).T
        rna2[design] = np.vstack(rna2[design]).T
        dna[design] = np.vstack(dna[design]).T
        assert rna1[design].shape == rna2[design].shape == dna[design].shape
    return names, seqs, rna1, rna2, dna, header

############################################################################
####### Write to files. ####################################################
    
def write_data(names, rna1, rna2, dna):
    header = ["ID", "DNA",
              "HepG2-RNA1", "K562-RNA1",
              "HepG2-RNA2", "K562-RNA2"]
    for design in designs:
        for p, promoter in enumerate(promoters):
            counts = np.hstack([names[design].reshape(-1, 1),
                                dna[design][:, p:p+1],
                                rna1[design][:, p*2:(p+1)*2],
                                rna2[design][:, p*2:(p+1)*2]]).astype('str')
            counts = np.vstack([header, counts])
            fn = "{}/{}-{}.counts".format(OUT_DIR, design, promoter)
            np.savetxt(fn, counts, delimiter = '\t', fmt = '%s')

if __name__ == '__main__':
    names, seqs, rna1, rna2, dna, header = read_data()
    write_data(names, rna1, rna2, dna)