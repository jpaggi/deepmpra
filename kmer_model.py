from sklearn.preprocessing import scale
from sklearn.model_selection import train_test_split
from sklearn.linear_model import SGDRegressor
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

def get_seqs(fn):
	seqs = {}
	with open(fn) as fp:
	    name = fp.readline().strip()[1:]
	    seq = ''
	    for line in fp:
	        if line[0] == '>':
	            seqs[name] = seq
	            name = line.strip()[1:]
	            seq = ''
	        else:
	            seq += line.strip().upper()
	seqs[name] = seq
	return seqs

def read_activities(fn):
	labels = {}
	with open(fn) as fp:
		fp.readline()
		for line in fp:
			name, act = line.strip().split()
			act = sum(map(float, act.split(':'))) / 2.0
			labels[name] = act
	return labels

##################################################################

BASES = ['A', 'C', 'G', 'T']

def kmerize(X, k):
    bases = ['00', '01', '10', '11']
    counts = []
    for seq in X:
        binary_seq = ''.join(map(lambda x: bases[x], seq))
        count = np.zeros((4**k,), dtype = int)
        for i in range(0, len(seq) - k + 1):
            count[int(binary_seq[i*2:(i+k)*2], 2)] += 1
        counts += [count]
    return np.vstack(counts)

def seqs_to_matrix(seqs):
    return np.vstack([map(lambda x: BASES.index(x), seq)
                     for seq in seqs])

def get_valid(labels, seqs, k):
    sequences = []
    y = []
    for i, (act_id, act) in enumerate(labels.items()):
        if act == float('inf'): continue
        y += [act]
        sequences += [seqs[act_id]]
    sequences = seqs_to_matrix(sequences)
    return kmerize(sequences, k), np.array(y)

#################################################################

def kmer_model(fasta, activities):
	seqs = get_seqs(fasta)
	activities = read_activities(activities)
	print len(seqs), len(activities)
	X, y = get_valid(activities, seqs, 6)
	y = scale(y)
	print X.shape, y.shape
	X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                   test_size = 0.2,
                                                   random_state = 42)
	model = SGDRegressor(penalty = 'l1').fit(X_train, y_train)
	pred = model.predict(X_test)
	print stats.spearmanr(y_test, pred)
	print stats.pearsonr(y_test, pred)
	plt.scatter(y_test, pred)
	plt.show()
	return model

##################################################################

def predict(seqs, model, k=6):
    seqs = kmerize(seqs_to_matrix(seqs), k)
    return model.predict(seqs)

def score_batch(batch, predict, genome, length):
	# Construct seqs
    ref_seqs, alt_seqs = [], []
    for chrom, pos, ref, alt, check in batch:
        ref_seqs += [genome.swapNs(
            genome.get_snv_seq(chrom, pos, ref, length, check))]
        alt_seqs += [genome.swapNs(
            genome.get_snv_seq(chrom, pos, alt, length, False))]
    # Run model
    ref_preds = predict(ref_seqs)
    alt_preds = predict(alt_seqs)
    return ref_preds, alt_preds


if __name__ == '__main__':
	import sys
	kmer_model(*sys.argv[1:])