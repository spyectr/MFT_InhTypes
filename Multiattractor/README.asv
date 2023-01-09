# Spiking networks with multiple inhibitory cell types

This repo runs simulation and mean field theory (MFT) analyses of recurrent networks of current-based leaky-integrate-and-fire neurons with multiple inhibitory cell types.
- Simulations are run with the scripts 'LIF_clustered_*.m'
- MFT runs with the 'MFT_*' scripts.

## MFT

### Attractor with 1 active cluster and the rest inactive

The simplest MFT calculation runs in 'MFT_selective.m' with the option options_main.a='sel' and calculates the Amit-Brunel type bifurcation plot
![bifurcation plot](https://github.com/spyectr/MFT_InhTypes/blob/master/Multiattractor/DATA/Results/Marcel2/MFT/figs/Marcel2_Pyr[g0]_FixedPoints_[6]_[14]_sel "Logo Title Text 1")

This script calculates the value of the firing rates as a function of the intra-cluster potentiation parameter J+. 
Below the critical value, all E clusters have the same firing rate. Above the critical value, the script computes the firing rate of 
the attractor with one cluster active at a high rate, and all other clusters inactive at a low rate (all of them grouped in one inactive 
population with the same firing rate).

### Map all multi-stable attractor 

The next MFT calculation runs in 'MFT_multistable.m' and starts by calculating the Amit-Brunel type bifurcation plot as above in a network with Q cluster, 
as a function of 
the intra-cluster potentiation parameter J+. The second MFT run uses the option options_main.a='all' and it searches for all the other stable attractors
with p clusters at a high rate and Q-p clusters at a low rate.


### Modulation of potential wells by external input

This script 'MFT_energylandscape.m' calculates the potential wells in control condition and with external inputs and plots them on top of each other.
The mode cue_mode=1 introduces the perturbation, whose parameters are in the structure 'cue', with fields:
- 'cue_option'={'Pyr'}; for perturbations targeting pyramidal cells. Other options: 'VIP', 'SST', 'PV'.
- 'cue_stat'={'gaussian'}; for perturbation modulating the quenched variance of the input to the chosen population. Other options: 'mean', for modulating the mean; 'noise' for introducing a dynamical input noise with mean 0 and specified variance, targeting the chosen population. 
- cueset={[0,0.1]}; For each entries in 'cue_option', specify the various value of the perturbation to be computed. This quantity is later set to the field 'cue.cue_value'.

Another option is nbins_grid=[80 80]; which is the size of the grid where the 2D energy landscape is calculated. Higher values mean finer grid, which leads to longer compute times.

The script calls the file 'MFT_run.m' which runs the following routines:
1. It first runs a MFT scan of the selective attractor with 1 cluster active as a function of J+, as in the above part.
2. It runs a second scan of MFT, calculating all the multistable attractors with multiple clusters active.

Steps 1. and 2. are computed for each value of the perturbation in 'cueset'.
3. Runs the potential well calculation for each value of 'cueset'.

Note: this script also evaluates the network linear response in the different attractor states at two values of J+: 'Jplus' and 'Jlow'.

Best practice is to first run the simple 'MFT_selective' script to find the pitchform bifurcation as a function of J+. Once the bifurcation is identified, set your parameters Jzero and Jstop to the right of the bifurcation and run a finer scan. Then run the potential well script.