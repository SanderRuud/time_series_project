#matplotlib.use('Agg')
import scipy.io
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.cluster import AgglomerativeClustering
from sklearn import cluster
import random

def create_plots(trans_matrix_csv, posterior_csv):
    matrix = pd.read_csv(trans_matrix_csv, index_col=None)
    transition_matrix = matrix[["V"+str(i) for i in range(1, len(matrix)+1)]].to_numpy()
    posterior = pd.read_csv(posterior_csv)
    fname = '../Mouse28-140313_BS0150_HMMready.mat'
    mat = scipy.io.loadmat(fname)
    angdata = np.ravel(np.array(mat['resampledAwakeHeadAngleData']))
    celldata = np.array(mat['celldata'])
    celldata = celldata.astype(int)
    ## toss out neurons that are not very active (because they boring)
    thrfrs = np.sum(celldata, 1)
    THR = 100 # I think this is a good threshold
    celldata = celldata[thrfrs>THR,:]
    Tfit = int(round(0.8*len(celldata[0,:])))
    Sfit = np.transpose(celldata[:,:Tfit])
    Stest = np.transpose(celldata[:,Tfit:])
    Sfit = Sfit.copy(order='C')
    Stest = Stest.copy(order='C')
    angfit = angdata[:Tfit]
    angtest = angdata[Tfit:]
    Sfit = (posterior[["S"+str(i) for i in range(1, len(matrix))]].to_numpy())
    state_sequence = posterior["state"].values
    used_states = posterior["state"].unique()
    k_means = cluster.KMeans(n_clusters=10, algorithm = "full")
    k_means.fit(transition_matrix)
    values = k_means.cluster_centers_
    f = defaultdict(list)

    for i, elem in enumerate(state_sequence):
        if i > 1200: break
        if not np.isnan(angfit[i]):
            elem = k_means.labels_[elem-1]
            f[elem].append(angfit[i])
    number_of_colors = len(f.keys())
    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                for i in range(number_of_colors+100)]
    for i, elem in enumerate(state_sequence):
        if i > 1200: break
        if not np.isnan(angfit[i]):
            elem = k_means.labels_[elem-1]
            plt.plot(i, angfit[i], ".", color=color[elem])
    plt.xlabel("time step")
    plt.ylabel("angle")
    plt.show()
    i=0
    for key, value in f.items(): 
        value = np.array(value)
        value = value[~np.isnan(value)]
        x = np.cos(value)
        y = np.sin(value)
        plt.scatter(np.average(x), np.average(y), label = key, s = 300, c=color[i]) # Plot average
        plt.scatter(x, y, label = key, c=color[i], s=.5) # Plot each point
        i+=1
    plt.show()

if __name__ == "__main__":
    create_plots("trans_matrix.csv", "posterior.csv") # All considered neurons
    create_plots("trans_matrix2.csv", "posterio2r.csv") # Only active neurons
    