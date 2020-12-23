import numpy as np
import scipy
from scipy import cluster

print("Quantizing")

def compute_thresholds(data, Q_BASE):
    thresholds = np.empty((sums.shape[2], Q_BASE-1), dtype=float)
    n_samples, n_time, n_genes = data.shape

    for i in range(n_genes):
        print("\r{0}%".format(np.around((i*100)/n_genes),4), end="")
        all_samples = data[:,:,i].flatten()
        codebook,distortion = scipy.cluster.vq.kmeans(all_samples, Q_BASE)
        codebook = np.sort(codebook)
        for q in range(Q_BASE-1):
            thresholds[i,q] = (codebook[q] + codebook[q+1])/2
    return thresholds



Q_BASE = 2

genedata = np.load(open("sorted_genes.npy","rb"))
genenames = np.load(open("sorted_names.npy","rb"))

sums = np.sum(genedata, axis=1)
thresholds = compute_thresholds(sums, Q_BASE)

quantized = np.ones(sums.shape, dtype=bool) * (Q_BASE-1)

for gene_i in range(sums.shape[2]):
    
    t = thresholds[gene_i]
    s = sums[:,:,gene_i]
    for q in range(Q_BASE-1):
        quantized[s < t[q], gene_i] = q

np.save(open("quantized_sorted_genes.npy", "wb"), quantized)
