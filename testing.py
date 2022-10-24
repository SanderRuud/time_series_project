
import sys
from os.path import dirname



#modulepath = "/Users/sanderruud/Downloads/temptidsrekker/pyhsmm-master"
modulepaths = ["/Users/sanderruud/Downloads/temptidsrekker/pyhsmm_spiketrains-master", "/Users/sanderruud/Downloads/temptidsrekker/pyhsmm_spiketrains-master/build",
            "/Users/sanderruud/Downloads/temptidsrekker/pyhsmm_spiketrains-master/pyhsmm_spiketrains",
            "/Users/sanderruud/Downloads/temptidsrekker/pyhsmm_spiketrains-master/pyhsmm_spiketrains/internals"]

#modulepaths = ["/Users/sanderruud/Downloads/temptidsrekker/pyhsmm_spiketrains-master"]
for module in modulepaths:
    sys.path.append(module)



print("Hello World!")
import pyhsmm
import pyhsmm_spiketrains.models


python -m "import pyhsmm_spiketrains.models"