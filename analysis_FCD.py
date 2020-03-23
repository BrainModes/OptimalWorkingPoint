# load necessary python modules
import scipy.io as sio
import numpy as np
import rpy2.robjects as rp
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

# define a function to calculate FCD
def FCD(BOLD, window_length, overlap):
    """
    BOLD singal   ....[regions * samples]
    window length ....[samples]
    overlap       ....[samples] 

    output: 
    FCD           ....[number_windows * number_windows]
    """
    import numpy as np
    n_regions = BOLD.shape[0]
    window_steps_size = window_length - overlap;
    n_windows = int(np.round((BOLD.shape[1] - window_length) / window_steps_size + 1))

    #compute FC for each window
    FC_t = np.zeros((n_regions,n_regions,n_windows))
    for i in range(n_windows):
        FC_t[:,:,i] = np.corrcoef(BOLD[:,window_steps_size*i:window_length+window_steps_size*i])
    
    
    # transform FC matrix into vector by just taking the upper triangle
    a,b =np.triu_indices(n_regions,1)
    tFC = FC_t[a,b,:].T
    
    
    # compute FCD, correlate the FCs with each other
    FCD = np.corrcoef(tFC);
    
    return FCD, FC_t

# calculate all empirical FCDs
subj_list = ['']          # specify a list of subject IDs here
path_to_empirical_fmri = ""  # path to find empirical fmri data
all_FCDs = np.zeros((len(subj_list),64,64)) # initialise array to store FCDs

# parameters for the FCD
window_length = 30; #[samples]
overlap       = 20; #[samples] 
for i in range(len(subj_list)):    
    #load empirical BOLD data from matlab format
    BOLD = sio.loadmat(path_to_empirical_fmri + '/' + subj_list[i] + '.mat')    
    #compute empirical FCD
    
    all_FCDs[i,:,:], _      = FCD(BOLD,window_length,overlap)
    #print(i)

# only 4 subjects were simulated with 22 min
# calculate simulated FCDs for each global coupling scaling factor
CSF = np.arange(0.025,0.0401,0.0001)
window_length = 30; #[samples]
overlap       = 20; #[samples]
all_simFCDs = np.zeros((4,151,63,63))
path_to_22min_simulation = ""
for i in range(len(subj_list)):
    for n in range(len(CSF)):
        mat = sio.loadmat(path_to_22min_simulation+subj_list[i]+"/"
                          +subj_list[ind[i]]+"_"+str(np.round(CSF[n],4))+".mat",
                         variable_names="Bold_data")
        BOLD = np.squeeze(mat['Bold_data'])
        BOLD = BOLD[40::4,:].T # cut out first 20s due to gradient, downsample to 0.5 Hz
        
        #compute sim FCD
        all_simFCDs[i,n,:,:], _      = FCD(BOLD,window_length,overlap)

# compare FCDs using the Kolmogorov-Smirnof distance measure
# import test from R
ks_dist = rp.r('ks.test')
all_KS_dist = np.zeros((4,151))
a,b = np.triu_indices(63,1)
for n in range(4):
    for i in range(151):
        all_KS_dist[n,i] = ks_dist(all_simFCDs[i,a,b],all_FCDs[i,a,b])[0][0]

# save results 
path_to_save = "" # specify where to save results
sio.savemat(path_to_save+"/FCD_KS_distance.mat",mdict={"all_KS_dist":all_KS_dist})