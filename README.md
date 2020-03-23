# OptimalWorkingPoint
This github repository contains simulation and analysis scripts for the paper 
"Identifying optimal working points of individual Virtual Brains: A large-scale brain
network modelling study" by Triebkorn et al. 


## Abstract 
Using The Virtual Brain (TVB, thevirtualbrian.org) simulation platform, we explored for 50 individual adult human
brains (ages 18-80), how personalized connectome based brain network modelling captures various empirical
observations as measured by functional magnetic resonance imaging (fMRI) and electroencephalography (EEG). We
compare simulated activity based on individual structural connectomes (SC) inferred from diffusion weighted imaging
with fMRI and EEG in the resting state. We systematically explore the role of the following model parameters:
conduction velocity, global coupling and graph theoretical features of individual SC. First, a subspace of the parameter
space is identified for each subject that results in realistic brain activity, i.e. reproducing the following prominent
features of empirical EEG-fMRI activity: topology of resting-state fMRI functional connectivity (FC), functional
connectivity dynamics (FCD), electrophysiological oscillations in the delta (3-4 Hz) and alpha (8-12 Hz) frequency
range and their bimodality, i.e. low and high energy modes. Interestingly, FCD fit, bimodality and static FC fit are
highly correlated. They all show their optimum in the same range of global coupling. In other words, only when our
local model is in a bistable regime we are able to generate switching of modes in our global network. Second, our
simulations reveal the explicit network mechanisms that lead to electrophysiological oscillations, their bimodal
behaviour and inter-regional differences. Third, we discuss biological interpretability of the Stefanescu-Jirsa-
Hindmarsh-Rose-3D model when embedded inside the large-scale brain network and mechanisms underlying the
emergence of bimodality of the neural signal.
With the present study, we set the cornerstone for a systematic catalogue of spatiotemporal brain activity regimes
generated with the connectome-based brain simulation platform The Virtual Brain.

## Content
### simulation script
TVB is a software for large scale brain network simulations. 
It is written in Python and can be downloaded for free from thevirtualbrian.org.
The following script makes use of this toolbox to run individual brain network simulations. 
The scripts take as an input the indivudual connectivity of a subject, which due to data protection laws cannot be shared openly. For empirical data and the full set of simulated data please contact Petra Ritter (petra.ritter@charite.de).
Two series of simulations were conducted in this study. One resulted in neural oscillations in the delta frequency range, the other in the alpha frequency range.
Comment or uncomment the indicated block in the script "simulation.py" to perform the different simulations. 
In the study we explore brain dynamics for 50 different subjects across a wide range of values for parameters global coupling and conduction speed. The present script only performs a single simulation for one combination of parameters. To explore the whole parameter space we used high performance computers and submitted multiple simulations with different parameters in parallel. 

### analysis scripts
The following scripts run in Python and use some R functions. 
Required toolboxes and modules: <br>
Python
1. numpy
2. scipy
3. rpy2 (to call R functions from Python)

R
1. stats (for the ks.test. i.e. Kolmogorov-Smirnoff test)
2. diptest (for Hartigan's diptest)

The "analysis.py" script, loops through the simulated time series and computes the fit between empirical and simulated functional connectivity, the dominant frequency and bimodality of the neural signal.
The "analysis_FCD.py" script, loops through the simulated time series and computes the fit between empirical and simulated functional connectivity using the Kolmogorov-Smirnoff distance.
