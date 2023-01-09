All network parameters can be found in the manuscript. Slight changes between the network with PV plasticity only
and PV & SST plasticity are denoted by an asterisk.

Please note that our background input are Poisson spike trains that elicit current changes in the receiving neuron.
Assemblies were trained by elevating the rate of Poisson spike trains that they received. The training phase ends after 370_000 ms.

Weight evolution of the network with PV and SST plasticity - weights_disinh.h5
Weight evolution of the network with PV plasticity only - weights_inh.h5

Nine averages weights were computed.

For exc synapses:
1 - intraassembly weights
2 - interassembly weights
3 - assemblies to background weights
4 - backgrounds to assembly weights
5 - background to background weights

For inh synapses:
6 - PV to assembly weights
7 - SST to assembly weights
8 - PV to background weights
9 - SST to background weights

3-5 are probably not important, and 6/8 and 7/9 could potentially be lumped together.

To see what I included in the files, simply write in Python:

import h5py
import numpy as np
with h5py.File(filename, 'r') as f:
    # List all groups
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]
    # Get the data
    data = list(f[a_group_key])

This will give you the names I saved the arrays under.

To access the arrays just write something like:

File = h5py.File('filename', 'r')
parameter = np.array(File.get('Name of array you want to access'))


