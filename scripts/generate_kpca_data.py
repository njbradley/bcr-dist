import sklearn.decomposition as skdecomp
import sklearn.manifold as skmanifold
from cmodule.build import bcrdist
import matplotlib.pyplot as plot
import numpy as np
import sys

array = bcrdist.dsbcellarray("../data/Hutchinson-BCR")

array.load_bd_data("../data/Hutchinson-BCRigh_DominantCDR3.csv", "../data/Hutchinson-BCRigl_DominantCDR3.csv")

array.generate_dist_matrix()
dist, ids = array.dist_matrix()

pcafunc = skdecomp.KernelPCA(75, kernel='precomputed')
tsnefunc = skmanifold.TSNE(n_components = 2, metric = 'precomputed')
tsne1dfunc = skmanifold.TSNE(n_components = 1, metric = 'precomputed')
normalized_dist = (dist / dist.max())
pcs = pcafunc.fit_transform(dist);
tsne = tsnefunc.fit_transform(dist);
tsne1d = tsne1dfunc.fit_transform(dist);


plot.scatter(pcs[:,0], pcs[:,1], s=1)
plot.savefig("../data/Hutchinson-BCR-pca-plot.png")
plot.clf()
plot.scatter(tsne[:,0], tsne[:,1], s=1)
plot.savefig("../data/Hutchinson-BCR-tsne-plot.png")

pcstxt = np.concatenate((np.array(ids).reshape(-1,1), pcs.astype(str)), axis=1)
np.savetxt("../data/Hutchinson-BCR-pcs.csv", pcstxt, delimiter=',', fmt='%s', comments='', header = "cell_index," + ','.join(["pc" + str(i) for i in range(75)]))

tsnetxt = np.concatenate((np.array(ids).reshape(-1,1), tsne.astype(str), tsne1d.astype(str)), axis=1)
np.savetxt("../data/Hutchinson-BCR-tsne.csv", tsnetxt, delimiter=',', fmt='%s', comments='', header = "cell_index,pre_tSNEx,pre_tSNEy,bcr_metric_tSNE")

alltxt = np.concatenate((np.array(ids).reshape(-1,1), pcs.astype(str), tsne.astype(str), tsne1d.astype(str)), axis=1)
np.savetxt("../data/Hutchinson-BCR-all-pcs-tsne-w1d.csv", alltxt, delimiter=',', fmt='%s', comments='', header = "cell_index," + ','.join(["pc" + str(i) for i in range(75)]) + ",pre_tSNEx,pre_tSNEy,bcr_metric_tSNE")
