# load necessary python modules
import numpy
import scipy
import os
# load TVB library
from tvb.simulator.lab import *
from tvb.basic.readers import FileReader


# set up the simulation

# select the SJHMR3D model as the local model and define its parameters
# for delta series choose 
oscilator = models.ReducedSetHindmarshRose(r=[0.006], a=[1.0], b=[3.0], c=[1.0], d=[5.0], s=[4.0], xo=[-1.6], K11=[0.5], K12=[0.1], K21=[0.15], sigma=[0.3], mu=[2.2], variables_of_interest=["xi","alpha"])
oscilator.configure()

# for alpha series uncomment the below code
# oscilator = models.ReducedSetHindmarshRose(r=[0.006], a=[1.0], b=[3.0], c=[1.0], d=[5.0], s=[4.0], xo=[-1.6], K11=[4.0], K12=[1.6], K21=[0.15], sigma=[0.4], mu=[2.2], variables_of_interest=["xi","alpha"])
# oscilator.configure()

# set up the structural connectivity
mypath=""  # insert the path to where all subjects connectivities are stored
subject="" # specify which subject to load

# load connectivity weights
reader = FileReader(file_path= mypath +'/' + subject +'/weights.txt')
w= reader.read_array(dtype=numpy.float64, skip_rows=0, use_cols=None, matlab_data_name=None)

# load region centers
reader = FileReader(file_path=mypath +'/' + subject +'/centres.txt')    
rl= reader.read_array( dtype="string", skip_rows=0, use_cols=(0,), matlab_data_name=None)
c= reader.read_array(dtype=numpy.float64, skip_rows=0, use_cols=(1, 2, 3), matlab_data_name=None)

# load connectivity tract lengths
reader = FileReader(file_path=mypath +'/' + subject +'/tract_lengths.txt')
tl= reader.read_array(dtype=numpy.float64, skip_rows=0, use_cols=None, matlab_data_name=None)

# confige the connectivity
white_matter = connectivity.Connectivity(region_labels=rl, weights=w, tract_lengths=tl, centres=c)        

# specify the coupling function
# here a linear scaling function was used
# csf --> global coupling scaling factor
# is one of the parameters explored in this study
# for delta series it was varied from 0.05 to 0.25 in steps of 0.01
# for alpha series it was varied from 0.025 to 0.04 in steps of 0.001
# in this simulation one one value is chosen, but loop across all values to reproduce the findings of the paper
csf = 0.05
white_matter_coupling = coupling.Linear(a=csf, b=[0.0])
white_matter_coupling.configure()

# specify a conduction speed to calculate time delays in the network 
# the 2nd parameter that was explored in this study
# for delta series it was varied from 20 to 100 in steps of 100
# for alpha series it was varied from 10 to 100 in steps of 10
# in this simulation one one value is chosen, but loop across all values to reproduce the findings of the paper
speed=100

# set up the integration scheme to solve the differential equations
# for delta series choose
heunint = integrators.HeunStochastic(dt=0.01220703125, noise=noise.Additive(nsig=array([1.0])))

# for alpha series uncomment the below code
heunint = integrators.HeunStochastic(dt=0.05, noise=noise.Additive(nsig=array([0.001])))


# specify what data to record from the simulation using tvb monitory
# choose subsampling of the neural signal
p=3.90625 #<-- 256Hz
momo = monitors.SubSample(period=p)

# choose to generate BOLD signal using a hemodynamic response function
hrffunction=equations.MixtureOfGammas()
pb=500
mama = monitors.Bold(period=pb, hrf_kernel=hrffunction)

# put both monitors together
what_to_watch = (momo, mama)

# configure the simulation
sim = simulator.Simulator(model=oscilator, connectivity=white_matter,
                          coupling=white_matter_coupling, conduction_speed=speed,
                          integrator=heunint, monitors=what_to_watch)

sim.configure()

# specify simulation length
# for delta series choose
sim_length = 180000

# for alpha series uncomment the below code
# sim_length = 300000

# for alpha series 22min simulation length use
sim_length = 1320000

# perform the simulation
subs_data = []
subs_time = []
bold_data = []
bold_time = []

for subs, bold in sim(simulation_length=sim_length):
    if not subs is None:
        subs_time.append(subs[0])
        subs_data.append(subs[1])
    
    if not bold is None:
        bold_time.append(bold[0])
        bold_data.append(bold[1])


# save the simulated time series
# declare a path where data should be saved
save_path=""
SUBS  = numpy.sum(numpy.array(subs_data),axis=3)
TSUBS = numpy.array(subs_time)
BOLD  = numpy.array(bold_data)
TBOLD = numpy.array(bold_time)
subsfile =save_path+'/'+subject+'csf_'+str(csf) +'_speed_' +str(speed)
scipy.io.savemat(subsfile, mdict={'subs_data': SUBS, 'subs_time': TSUBS ,'bold_data': BOLD, 'bold_time': TBOLD })





