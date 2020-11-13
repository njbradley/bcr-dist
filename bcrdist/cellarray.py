import sklearn.decomposition as skdecomp
import sklearn.manifold as skmanifold
import sklearn.preprocessing as skpre
from . import cbcrdist
import matplotlib.pyplot as plot
import numpy as np
import sys
import os

cbcrdist.init(os.path.dirname(os.path.abspath(__file__)) + '/')

class cellarray(cbcrdist.bcellarray):
    '''
    An array of bcells, with either a light and heavy chain or a single heavy chain
    This python subclass of the pure c class bcrdist.dsbcellarray implements
    new plotting functions like generate kpca data, which invoke python
    libraries and functions
    '''
    
    def tsneplot(self, path = None, colorby = None):
        if (path == None):
            path = self.name() + "-tsneplot.png"
        
        tsnefunc = skmanifold.TSNE(n_components = 2, metric = 'precomputed')
        
        dist, ids = self.distmatrix()
        data = self.tolist()
        
        tsne = tsnefunc.fit_transform(dist)
        
        
        if (type(colorby) == type(lambda x: x)):
            colors = [colorby(cell) for cell in data]
        elif (colorby == "clonotype"):
            clonotypes = [cell[1] for cell in data]
            print (clonotypes)
            clone_indexes = list(set(clonotypes))
            clonotypes.insert(0, 'None')
            colors = [clone_indexes.index(i) for i in clonotypes]
            print (colors)
        
        
        plot.clf()
        plot.title(self.name() + " tSNE")
        if (colorby == None):
            plot.scatter(tsne[:,0], tsne[:,1], s=1)
        else:
            plot.scatter(tsne[:,0], tsne[:,1], s=1, c=colors)
            plot.colorbar()
        plot.savefig(path)
    
    def savePCs(self, path = None):
        if (path == None):
            path = self.name() + "-pcs.csv"
        dist, ids = self.distmatrix()
        pcafunc = skdecomp.KernelPCA(75, kernel='precomputed')
        kernel = np.exp(-dist**2 / dist.max()**2)
        pcs = pcafunc.fit_transform(kernel)
        
        pcstxt = np.concatenate((np.array(ids).reshape(-1,1), pcs.astype(str)), axis=1)
        np.savetxt(path, pcstxt, delimiter=',', fmt='%s', comments='', header = "cell_index," + ','.join(["pc" + str(i) for i in range(75)]))
    
    def umapplot(self, path=None, colorby=None):
        if (path == None):
            path = self.name() + "-umapplot.png"
        
        import umap
        dist, ids = self.distmatrix()
        data = self.tolist()
        
        umapfunc = umap.UMAP(metric = "precomputed")
        umapout = umapfunc.fit_transform(dist)
        
        if (type(colorby) == type(lambda x: x)):
            colors = [colorby(cell) for cell in data]
        elif (colorby == "clonotype"):
            clonotypes = [cell[1] for cell in data]
            counts = [(clonotypes.count(i),i) for i in list(set(clonotypes))]
            counts.sort(reverse=True)
            print(counts[:20])
            
            clone_indexes = [i[1] for i in counts[:20]]
            colors = [(-1 if not i in clone_indexes else clone_indexes.index(i)) for i in clonotypes]
            
        
        plot.clf()
        plot.title(self.name() + " UMAP")
        if (colorby == None):
            plot.scatter(umapout[:,0], umapout[:,1], s=1)
        else:
            plot.scatter(umapout[:,0], umapout[:,1], s=1, c=colors)
            plot.colorbar()
        plot.savefig(path)
    
    def generate_kpca_data(self):
        dist, ids = self.distmatrix()
        
        pcafunc = skdecomp.KernelPCA(75, kernel='precomputed')
        tsnefunc = skmanifold.TSNE(n_components = 2, metric = 'precomputed')
        tsne1dfunc = skmanifold.TSNE(n_components = 1, metric = 'precomputed')
        
        #kernel = 1 - (dist / dist.max())
        kernel = np.exp(-dist**2 / dist.max()**2)
        pcs = pcafunc.fit_transform(kernel)
        tsne = tsnefunc.fit_transform(dist)
        tsne1d = tsne1dfunc.fit_transform(dist)
        #umapin = skpre.StandardScaler().fit_transform(dist)
        
        
        plot.scatter(pcs[:,0], pcs[:,1], s=1)
        plot.savefig(self.name() + "-pca-plot.png")
        plot.clf()
        plot.scatter(tsne[:,0], tsne[:,1], s=1)
        plot.savefig(self.name() + "-tsne-plot.png")
        
        # umapfunc = umap.UMAP()
        # umapout = umapfunc.fit_transform(skpre.StandardScaler().fit_transform(pcs))
        # tsnefunc = skmanifold.TSNE(n_components = 2)
        # tsne = tsnefunc.fit_transform(dist)
        #
        # plot.clf()
        # plot.scatter(umapout[:,0], umapout[:,1], s=1)
        # plot.savefig(self.name() + "-umap-plot-pcs.png")
        # plot.clf()
        # plot.scatter(tsne[:,0], tsne[:,1], s=1)
        # plot.savefig(self.name() + "-tsne-plot-pcs.png")
        
        pcstxt = np.concatenate((np.array(ids).reshape(-1,1), pcs.astype(str)), axis=1)
        np.savetxt(self.name() + "-pcs.csv", pcstxt, delimiter=',', fmt='%s', comments='', header = "cell_index," + ','.join(["pc" + str(i) for i in range(75)]))
        
        tsnetxt = np.concatenate((np.array(ids).reshape(-1,1), tsne.astype(str), tsne1d.astype(str)), axis=1)
        np.savetxt(self.name() + "-tsne.csv", tsnetxt, delimiter=',', fmt='%s', comments='', header = "cell_index,pre_tSNEx,pre_tSNEy,bcr_metric_tSNE")
        
        alltxt = np.concatenate((np.array(ids).reshape(-1,1), pcs.astype(str), tsne.astype(str), tsne1d.astype(str)), axis=1)
        np.savetxt(self.name() + "-all-pcs-tsne-w1d.csv", alltxt, delimiter=',', fmt='%s', comments='', header = "cell_index," + ','.join(["pc" + str(i) for i in range(75)]) + ",pre_tSNEx,pre_tSNEy,bcr_metric_tSNE")

        print( "sucessfully saved tsne and kpca data to file " + self.name() + "-all-pcs-tsne-w1d.csv" )
    
    def __repr__(self):
        message = "<bcrdist.cellarray object\n"
        message += self.summary()
        message += ">"
        return message
