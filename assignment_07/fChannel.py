import numpy as np

# ......

# fChannel.py

# ........

def myf_channel(sys_param, param_channel):

    ##parameters##

    Nt = sys_param["Nt"]

    Nr = sys_param["Nr"]

    seed_seq = param_channel["seed_seq"]


    rng = np.random.default_rng(seed = seed_seq)

    
    channel = {}

    #Large-scale pathloss

    PL = 1

    ##Small-scale fading ###

    H_mat = rng.normal(loc = 0, scale = np.sqrt(1/2), size = (Nr, Nt)) \
                   + 1j * rng.normal(loc=0, scale= np.sqrt(1/2), size=(Nr,Nt))
    
    ## Final channel ####

    H_mat = np.sqrt(PL)  * H_mat


    ##SVD ###

    U_mat, S_array, V_mat_H = np.linalg.svd(H_mat, full_matrices=True)

    channel["H_mat"] = np.copy(H_mat)

    channel["U_mat"] = np.copy(U_mat)

    channel["S_mat"] = np.copy(np.diag(S_array))

    channel["V_mat"] = np.transpose(np.conjugate(V_mat_H))

    return channel

