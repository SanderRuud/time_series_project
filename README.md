# Uncovering Neural States in Rodents

This repository contains the code and data used in a semester project in TMA4285 - Time Series. The report can be read [here](report_time_series_project.pdf). The report will be graded as approved or not-approved.

Running the supplied code (runHMMonmousedata.py / shortpostprocess.py) is not necessary, but an installation guide is provided below.




# Installation of required packages 

The packages are tedious, some of them aren't updated.  


Adding some packages to a new virtual environment. 


```
conda create --name time_series python=3.7 -y
conda activate time_series
conda install -y pip numpy scipy matplotlib cython nose requests networkx
```
Either  clone the repos and, ```cd {path}```, ```pip install .``` and add to ```sys.path``` , or download directly using the github links.

```
pip install git+https://github.com/mattjj/pybasicbayes 
pip install git+https://github.com/mattjj/pyhsmm
python -c "import pyhsmm"
```
Depending on what compiler version for C++ you have, you might get an error for ```pip install git+https://github.com/slinderman/pyhsmm_spiketrains```

Specify the versions of CXX and CC, either with ```CXX=g++-12 CC=gcc-12 pip install git+https://github.com/slinderman/pyhsmm_spiketrains```, or clone the repo. 
```
CXX=g++-12 CC=gcc-12 python setup.py install build_ext --inplace --user
```
Then, install hips

```
pip install git+https://github.com/HIPS/hips-lib
```
Then fix some things that are no longer allowed in python 3.7

Note that ```pip install scipy==1.1.0``` avoids ```logsumexp``` being moved from ```scipy.misc``` to ```scipy.special```. (Alternatively, change the import in the hips-package)

**pyhsmm_spiketrains/internals/poisson_observations.py [53]**  (Sublist parameters are not supported in Python 3.x).

```
@hypers.setter # Old
def  hypers(self, (alpha_0, beta_0)):

@hypers.setter # New
def  hypers(self, v):
	alpha_0, beta_0 = v
```

**pyhsmm_spiketrains/internals/poisson_observations.py [78]** 

```
if  None  not  in (n,tots): # Old
	for  p, tot  in  zip(self.poissons,tots):

if (n  is  not  None) and (tots  is  not  None): # New
	for  p, tot  in  zip(self.poissons,tots):
```

**pyhsmm_spiketrains/models.py [306]**

```
for  n  in  xrange(N): # Old

for  n  in  range(N): # New
```
Then runHMMonmousedata.py can be run.
