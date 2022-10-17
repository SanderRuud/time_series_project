import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

from pyhsmm.basic.distributions import PoissonDuration
from pybasicbayes.util.text import progprint_xrange

import scipy.io

import pyhsmm_spiketrains.models
reload(pyhsmm_spiketrains.models)

import sys


# Generate a synthetic dataset
N_iter = 300    # Number of iterations of Gibbs sampling
## in the paper they needed a ton of iterations for it to work!!!!!

#fname = sys.argv[1]
## data binned into 150ms bins (this number is just pulled out of a hat)
fname = 'Mouse28-140313_BS0150_HMMready.mat'
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


basename = '%s_project_HMM_R%05d'%(fname[:(-4)], int(round(np.random.rand()*10000)))

N = len(Sfit[:,0])


alpha_obs = 1.0
beta_obs = 1.0
priorhyperparams = 5.0
hmmmodel = pyhsmm_spiketrains.models.PoissonHDPHMM(N=N, K_max=100, alpha_a_0=priorhyperparams, alpha_b_0=1.0, gamma_a_0=priorhyperparams, gamma_b_0=1.0, init_state_concentration=1.0, alpha_obs=alpha_obs, beta_obs=beta_obs)
hmmmodel.add_data(Sfit)

# Fit the test model with Gibbs sampling
loglikes_fit = []
loglikes_test = []
count = 0
for itr in progprint_xrange(N_iter):
    hmmmodel.resample_model()
    loglikes_fit.append(hmmmodel.log_likelihood(Sfit))
    loglikes_test.append(hmmmodel.log_likelihood(Stest))
    count += 1

    if(np.mod(count,50)==0):
      # Get the inferred state sequence
      hmmmodel.relabel_by_usage()
      Z_train_inf = hmmmodel.stateseqs[0]
      N_used_inf = 0
      used_states = []
      for i in hmmmodel.used_states:
        used_states.append(i)
        N_used_inf += 1
      print('Number of states used', N_used_inf)

      outy = {}
      outy['transition_matrix'] = hmmmodel.A[:N_used_inf, :N_used_inf]
      outy['used_states'] = used_states
      outy['state_sequence'] = hmmmodel.stateseqs[0]
      outy['Sfit'] = Sfit
      outy['angfit'] = angfit
      scipy.io.savemat('%s_post_hmm_ITER%04d.mat'%(basename, count), outy)



# Plot the log likelihood over time
plt.figure()
plt.plot(loglikes_fit, 'b')
plt.plot([0,N_iter], hmmmodel.log_likelihood(Sfit) * np.ones(2), ':k')
plt.xlabel("Iteration")
plt.ylabel("Log Likelihood")
plt.savefig('%s_log_like.pdf'%basename)

# Plot the log likelihood over time
plt.figure()
plt.plot(loglikes_test, 'r')
plt.xlabel("Iteration")
plt.ylabel("Log Likelihood")
plt.savefig('%s_log_like_test.pdf'%basename)

# Visualize the data and the inferred state sequences
plt.figure()
plt.subplot(211)
plt.imshow(Sfit.T[:,:100], interpolation="none", cmap="Greys", vmin=0, vmax=Sfit.max())
plt.title("Spike train")
plt.subplot(212)
plt.title("Inferred states")
plt.imshow(Z_train_inf.reshape((1,-1))[:,:100], aspect=10.0, cmap="Paired", interpolation="none", vmin=0, vmax=N_used_inf)
plt.savefig('%s_state_seqs.pdf'%basename)

# Visualize the inferred transition matrices
plt.figure()
plt.imshow(hmmmodel.A[:N_used_inf, :N_used_inf], interpolation="none", cmap="Greys", vmin=0, vmax=1)
plt.title("Inf. Transitions")
plt.savefig('%s_transition_matrix.pdf'%basename)

plt.show()
