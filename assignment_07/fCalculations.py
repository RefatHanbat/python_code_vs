import time

import numpy as np

# from scipy.linalg import block_diag

# ,,,

# fCalculations.py

# ,,,

def myf_constellation(sys_param):

    constellation_type = sys_param["constellation_type"]

    ##Function###

    E_symbol_avg = 1 # Unit average symbol energy

    constellation = {}

    if(constellation_type=="BPSK"):

        constellation_symbols = np.sqrt(E_symbol_avg/1)\
             *np.array([[-1],\
                       [+1]])
        constellation_bits = np.array([[0],\
                                       [1]])
    elif(constellation_type == "QPSK"):

        constellation_symbols = np.sqrt(E_symbol_avg/2)\
        *np.array([ [(-1) + 1j * (-1)],\
                   
                   [(-1) + 1j * (+1)],\
                   
                   [(+1) + 1j * (-1)],\
                   [(+1) + 1j * (+1)] ])
        
        constellation_bits = np.array([[0,0],\
                                       [0,1],\
                                        [1,0],\
                                        [1,1] ])
        
    constellation["constellation_symbols"] = constellation_symbols

    constellation["constellation_bits"] = constellation_bits

    return constellation


def myf_symbol_indices(sys_param, param_symbol_indices):

    Ns = sys_param["Ns"]

    constellation_size = sys_param["constellation_size"]

    seed_seq = param_symbol_indices["seed_seq"]

    ## Function###

    ##Random seed

    rng = np.random.default_rng(seed=seed_seq)

    symbol_indices = rng.integers(low=0, high=constellation_size, size=Ns)

    return symbol_indices

def myf_x_vec(sys_param, constellation,symbol_indices):

    ##parameters##

    Ns = sys_param["Ns"]

    constellation_symbols = constellation["constellation_symbols"]


    ###Function##

    x_vec = np.zeros((Ns,1), dtype= np.complex64)


    for ind in range(0,Ns):
        
        x_vec[ind][0] = constellation_symbols[symbol_indices[ind]][0]

        

    return x_vec


def myf_P(sys_param, channel, mode="EPA"):

    ##parameters##

    Nt = sys_param["Nt"]

    Ns = sys_param["Ns"]

    Es = sys_param["Es"]

    NO = sys_param["NO"]

    S_mat = channel["S_mat"]

    V_mat = channel["V_mat"]
    
    ##Function##

    if(mode=="EPA"):

        P_mat = np.sqrt(Es/Ns)*np.vstack((np.eye(Ns, dtype=np.float32),np.zeros((Nt-Ns, Ns),\
                                                                            dtype=np.float32)))

        P_mat_reduced = np.copy(P_mat)

        Ns_opt = np.copy(Ns)

    elif(mode=="SVD"):

        Q_mat_1_opt = np.zeros((Ns,Ns), dtype=np.float32)

        Ns_opt = 1

        C_opt = 0

        for ind1 in range(0,Ns):

            flag_some_q_neg = 0

            Q_mat_1_temp = np.zeros((Ns,Ns), dtype=np.float32)

            Ns_temp  = ind1 + 1

            C_temp = 0

            sum_temp = 0

            for ind2 in range(0, Ns_temp):

                sum_temp = sum_temp + NO / S_mat[ind2][ind2]



            mu_temp = (1/Ns_temp) * (Es + sum_temp)

            for ind2 in range(0,Ns_temp):

                Q_mat_1_temp [ind2][ind2] = mu_temp - NO/S_mat[ind2][ind2]

                if(Q_mat_1_temp[ind2][ind2] < 0 ):

                    flag_some_q_neg = 1

                    break
            if(flag_some_q_neg != 1 ):

                for ind2 in range(0,Ns_temp):

                    C_temp = C_temp + np.log2(1 + (1/NO)* (S_mat[ind2][ind2]**2)* Q_mat_1_temp[ind2][ind2])

                    if(C_temp > C_opt):

                        Q_mat_1_opt = np.copy(Q_mat_1_temp)

                        Ns_opt = np.copy(Ns_temp)

                        C_opt = np.copy(C_temp)

        V_mat_1 = V_mat[:, 0:Ns]

        P_mat = V_mat_1@np.sqrt(Q_mat_1_opt)

        Q_mat_1_opt_reduced = Q_mat_1_opt[:Ns_opt, :Ns_opt]

        V_mat_1_reduced = V_mat[:,0:Ns_opt]

        P_mat_reduced = V_mat_1_reduced @ np.sqrt(Q_mat_1_opt_reduced)
    
    return P_mat, P_mat_reduced,Ns_opt

def myf_y_vec(sys_param, channel, param_y_vec):

    Nr = sys_param["Nr"]

    NO = sys_param["NO"]

    H_mat = channel["H_mat"]

    seed_seq = param_y_vec["seed_seq"]

    s_vec = param_y_vec["s_vec"]

    ##Function##

    ##Random seed

    #np.random.RandomState(seed=0)

    rng = np.random.default_rng(seed = seed_seq)

    z_vec = rng.normal(loc=0, scale=np.sqrt(NO/2), size=(Nr,1))\
            +1j * rng.normal(loc= 0, scale= np.sqrt(NO/2), size=(Nr,1))
    
    y_vec = H_mat @s_vec + z_vec

    return y_vec

def myf_C_Pandw(sys_param, channel, P_mat, Ns_opt, solutions):

    Nr = sys_param["Nr"]

    NO = sys_param["NO"]\
    
    H_mat = channel["H_mat"]

    W_mat = solutions["W_mat"]

    ##Function ###

    R_z_mat = NO * np.eye(Nr,dtype=np.float32)

    P_mat_H = np.transpose(np.conjugate(P_mat))

    W_mat_H = np.transpose(np.conjugate(W_mat))

    H_mat_H = np.transpose(np.conjugate(H_mat))

    c = np.log2(np.linalg.det(np.eye(Ns_opt, dtype=np.float32)+\
                        np.linalg.inv(W_mat@R_z_mat@W_mat_H)@W_mat@H_mat@P_mat@P_mat_H@H_mat_H@W_mat_H))
    

    return c

def myf_Num_errors(sys_param, constellation, solutions, param_Num_errors):

    constellation_type = sys_param["constellation_type"]

    constellation_bits = constellation["constellation_bits"]

    W_mat = solutions["W_mat"]

    symbol_indices = param_Num_errors["symbol_indices"]

    Ns_opt = param_Num_errors["Ns_opt"]

    y_vec = param_Num_errors["y_vec"]


    ##Estimated x_vec

    x_vec_est = W_mat @ y_vec

    

    ##Estimated symbol indices

    symbol_indices_dec = np.zeros(Ns_opt)

    for ind in range(0,Ns_opt):
        
        if(constellation_type=="BPSK"):

            if(np.real(x_vec_est[ind][0]) < 0 ):

                symbol_indices_dec[ind] = 0

            elif(np.real(x_vec_est[ind][0]) >= 0 ):

                symbol_indices_dec[ind] = 1

        elif(constellation_type=="QPSK"):

            if((np.real(x_vec_est[ind]) < 0) and ((np.imag(x_vec_est[ind][0]) < 0))):

                symbol_indices_dec[ind] = 0

            elif((np.real(x_vec_est[ind]) < 0) and ((np.imag(x_vec_est[ind][0]) >= 0))):

                symbol_indices_dec[ind] = 1

            elif((np.real(x_vec_est[ind]) >= 0) and ((np.imag(x_vec_est[ind][0]) < 0))):

                symbol_indices_dec[ind] = 2

            elif((np.real(x_vec_est[ind]) >= 0) and ((np.imag(x_vec_est[ind][0]) >= 0))):

                symbol_indices_dec[ind] = 3

                
    Num_errors_symbol = 0

    Num_errors_bit = 0

    for ind in range(0,Ns_opt):

        if(symbol_indices[ind] != symbol_indices_dec[ind]):

            Num_errors_symbol = Num_errors_symbol + 1

            Num_errors_bit = Num_errors_bit + 1

            Num_errors_bit = Num_errors_bit + np.sum(np.abs(constellation_bits[int(symbol_indices[ind])]\
                            - constellation_bits[int(symbol_indices_dec[ind])]))
            
    Num_errors = {}

    Num_errors["Num_errors_symbol"] = Num_errors_symbol

    Num_errors["Num_errors_bit"] = Num_errors_bit
            
    return Num_errors
