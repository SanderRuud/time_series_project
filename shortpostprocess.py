import matplotlib
matplotlib.use('Agg')
import scipy.io
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
import pylab as pltt
import sys

import pickle

fname = sys.argv[1]
basename = fname[:-4]
mat = scipy.io.loadmat(fname)
#['used_states', 'state_sequence', 'Sfit', '__header__', 'transition_matrix', '__globals__', 'angfit', '__version__']

transition_matrix = np.array(mat['transition_matrix'])

angfit = np.ravel(np.array(mat['angfit']))
state_sequence = np.ravel(np.array(mat['state_sequence']))
Sfit = np.array(mat['Sfit'])

#print shape(state_sequence), shape(Sfit)

used_states = np.ravel(np.array(mat['used_states']))

meanvalues = np.zeros(len(used_states))
stdvalues = np.zeros(len(used_states))
counts = np.zeros(len(used_states))
numinstances = np.zeros(len(used_states))
tot = 0
for i in used_states:
  count = np.sum(state_sequence==i)
  tot += count

  angs = angfit[state_sequence==i]
  if(sum(np.isnan(angs))>0):
    vals = 1.*(state_sequence==i)
    vals = vals[1:]-vals[:(-1)]
    numinstances[i] = np.sum(vals>0.5) ## i think this is correct

    angs = angs[~np.isnan(angs)]
    if(sum(angs<0.)>0):
      angs[angs<0.] += 2.*pi

    #### Note only filtering here because these are the ones we can check
    if(len(angs)>10):
      meanvalues[i] = scipy.stats.circmean(angs)
      stdvalues[i] = scipy.stats.circstd(angs)
      counts[i] = count


print(fname)
inds = np.argsort(meanvalues)
meanvalues = meanvalues[inds]
stdvalues = stdvalues[inds]
counts = counts[inds]
numinstances = numinstances[inds]

## resort according to angle
sortedTM = np.zeros(np.shape(transition_matrix))
sortedTM[:] = np.nan
for i in range(len(inds)):
  if(counts[i] < 10):
    continue
  for j in range(len(inds)):
    if i==j:
      continue
    if(counts[j] < 10):
      continue
    sortedTM[i,j] = (transition_matrix[inds[i], inds[j]])

for i in range(len(inds)):
  if(counts[i]>-1):
    print(i+1, meanvalues[i], stdvalues[i], counts[i], numinstances[i])


## remove the empty ones
guys = []
for i in range(len(inds)):
  if(stdvalues[i] == 0):
    continue
  for j in range(len(inds)):
    if i==j:
      continue
    if(stdvalues[j] == 0):
      continue
    guys.append( (i+1, j+1, {'myweight':( (transition_matrix[inds[i], inds[j]]) )} ) )


import networkx as nx
G = nx.Graph()
G.add_edges_from(guys)
pos = nx.spectral_layout(G,weight='myweight') #,  iterations=10000)
nx.draw_networkx(G,pos)
pltt.savefig('%s__spectral_layout.png'%(basename))


#plt.figure()
#plt.imshow(sortedTM, cmap='jet', interpolation='nearest')
