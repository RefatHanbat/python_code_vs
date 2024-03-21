import time 

import numpy as np

def myf_algorithm_1(sys_param, channel, P_mat, Ns_opt, mode="RxMF"):

    Nr = sys_param["Nr"]

    NO = sys_param["NO"]

    H_mat = channel["H_mat"]

    ##Function ###

    R_z_mat = NO * np.eye(Nr,dtype=np.float32)

    R_z_mat_inv = np.linalg.inv(R_z_mat)

    P_mat_H = np.transpose(np.conjugate(P_mat))

    H_mat_H = np.transpose(np.conjugate(H_mat))

    if(mode == "RxMF"):

        W_mat = P_mat_H @ H_mat_H @ R_z_mat_inv

    elif(mode == "RxZF"):

        W_mat = np.linalg.inv(P_mat_H @ H_mat_H @R_z_mat_inv@H_mat@P_mat)@P_mat_H @H_mat_H@R_z_mat_inv

    elif(mode =="RxMMSE"):

        W_mat = np.linalg.inv(np.eye(Ns_opt, dtype=np.float32) + 
                P_mat_H@H_mat_H@R_z_mat_inv@H_mat@P_mat)@P_mat_H@H_mat_H@R_z_mat_inv
        
    solutions = {}

    solutions["W_mat"] = np.copy(W_mat)

    return solutions

def myf_algorithm_2(sys_param, channel, mode="svd"):

    Ns = sys_param["Ns"]

    U_mat = channel["U_mat"]

    
    if(mode == "SVD"):
        
        W_mat = np.transpose(np.conjugate(U_mat[:, 0:Ns]))

    solutions = {}

    solutions["W_mat"] = np.copy(W_mat)

    return solutions