# load necessary python modules
import numpy as np
import scipy as scp
import scipy.io
import scipy.signal
import rpy2.robjects as rp

subject = '';  # specify subject to analyse
speed = 10     # specify which speed to analyse, this was done to split analysis across multiple cores


# specify where all empirical functional connectivities can be found
path_to_all_empirical_FCs 
# load the empirical FC as a numpy array, i.e. a 2D matrix with dimensions [68,68]
empFC      = scipy.io.loadmat(path_to_all_empirical_FCs)

# initialis arrays to store analysis data
all_simFC  = np.zeros((151,68,68)) # to store simulated FCs
sim_empFC  = np.zeros((151,1))     # to store the results of simulated to empirical FC comparison

avg_bimod = np.zeros((151,4))     # store statistics of bimodality tests for the average across all region signal, csf*[ p_dip, stat_Dip, p_dip_log, stat_Dip_log]
avg_freq  = np.zeros((151))       # store dominant frequency for the average across all region signal
reg_bimod = np.zeros((151,68,4))  # store statistics of bimodality tests for each region, csf*region*[p_dip, stat_Dip, p_dip_log, stat_Dip_log]
reg_freq  = np.zeros((151,68))    # store dominant frequency for each region

# global coupling scaling factor to loop across
CSF = np.arange(0.025,0.0401,0.0001)

# importing Hartigan's diptest for bimodality from R
# this package needs to be installed in the R distribution
d   = rp.r('diptest::dip.test')

# specify where simulations are stored
sim_results_path = ''
for csf in range(len(CSF)):
    
    #load simulated timeseries, subsampled neural signal and BOLD/fMRI
    subs = scp.io.loadmat(sim_results_path + subject +'_csf_'+str(CSF[csf])+'_speed_'+str(speed)+'.mat')['subs_data'].mean(axis=1)
    subs = subs[2*200:,:]
    bold = scp.io.loadmat(sim_results_path + subject +'_csf_'+str(CSF[csf])+'_speed_'+str(speed)+'.mat')['Bold_data']
    bold = np.squeeze(bold[40:,:,:,:].mean(axis=(1,3))).T
    
    # calculate simulated to empirical FC fit
    all_simFC[csf,:,:] = np.corrcoef(bold)
    
    sim_empFC[csf] = np.corrcoef(all_simFC[csf,:,:].flatten(),empFC.flatten())[0,1]
    
    # calculate statistics for average signal
    f, t, Sxx  = scp.signal.spectrogram(subs.mean(axis=1), fs=200, window=('tukey', 0.25), nperseg=128, noverlap=110, nfft=4*200)
    mean_power = np.mean(Sxx, axis=1)
    # dominant frquency
    avg_freq[csf]   = f[np.argmax(mean_power)]

    # bimodality
    dip = d(rp.FloatVector((Sxx[np.argmax(mean_power),:])))
    avg_bimod[csf,0] = dip[1][0]
    avg_bimod[csf,1] = dip[0][0]
    
    # now for log(power)
    dip = d(rp.FloatVector((np.log(Sxx[np.argmax(mean_power),:]))))
    avg_bimod[csf,2] = dip[1][0]
    avg_bimod[csf,3] = dip[0][0]
    
    
    # calculate statistis for each region
    for i in range(68):
        # dominant frequency
        f, t, Sxx = scp.signal.spectrogram(subs[:,i], fs=200, window=('tukey', 0.25), nperseg=128, noverlap=110, nfft=4*200)
        mean_power = np.mean(Sxx, axis=1)
        reg_freq[csf,i]   = f[np.argmax(mean_power)]

        # bimdality test                       
        dip = d(rp.FloatVector((Sxx[np.argmax(mean_power),:])))
        reg_bimod[csf,i,0] = dip[1][0]
        reg_bimod[csf,i,1] = dip[0][0]
        
        #now for log(power)
        dip = d(rp.FloatVector((np.log(Sxx[np.argmax(mean_power),:]))))
        reg_bimod[csf,i,2] = dip[1][0]
        reg_bimod[csf,i,3] = dip[0][0]


#save results, specify a path here
file = save_path+subject+"_speed_"+str(speed)+".mat"
scipy.io.savemat(file,mdict={'all_simFC': all_simFC, 'sim_empFC':sim_empFC, 'avg_bimod':avg_bimod, 
                             'avg_freq': avg_freq, 'reg_bimod': reg_bimod, 'reg_freq ':reg_freq})
