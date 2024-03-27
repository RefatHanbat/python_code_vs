import time
import numpy as np
from fCalculations import *

# .......
# fAlgorithms.py
# .......

def myf_algorihtm_1(sys_param, channel):
    
    P_S = sys_param["P_S"]
    
    P_D_bar = sys_param["P_D_bar"]

    r_P_bar = sys_param["r_P_bar"] 

    No_W = sys_param["No_W"]

    No_W_unc = sys_param["No_W_unc"]

    No_D = sys_param["No_D"]

    No_R = sys_param["No_R"]

    err_min = sys_param["err_min"]

    h_SD = channel["h_SD"]

    h_SR = channel["h_SR"]

    h_DR = channel["h_DR"]

    h_DW = channel["h_DW"]

    h_DD_tilde = channel["h_DD_tilde"]

    if(r_P_bar == 0 ):

        P_D_ub = np.array([((No_W_unc**(1-4*err_min)-1/No_W_unc)*No_W)/(np.abs(h_DW)**2),P_D_bar])

        

    else:

        # P_D_ub = np.array([\
        #                  (1/np.abs(h_DR)**2) * ((np.abs(h_SR)**2 * P_S/(2**r_P_bar - 1))-No_R),\
        #                   (1/np.abs(h_DD_tilde)**2) * ((np.abs(h_SD)**2 * P_S/(2**r_P_bar - 1))-No_D),\
        #                   ((No_W_unc**(1-4 * err_min) - 1/(No_W_unc)) - No_W/(np.abs(h_DW)**2) ),\
        #                   P_D_bar
        #                             ])

        # P_D_ub = np.array([\
        #                 ((np.abs(h_SR)**2 * P_S) / (2**(r_P_bar) - 1) -
        #         No_R) (np.abs(h_DR)**2),\
        #                 ((np.abs(h_SD)**2 * P_S) / (2**(r_P_bar) - 1) -
        #         No_D) / (np.abs(h_DD_tilde)**2),\
        #                  ((No_W_unc**(1 - 4 * err_min) - 
        #         1 / No_W_unc) * No_W) / (np.abs(h_DW)**2),\
        #                                 P_D_bar
        #                                           ])
        h_SR_sq_mgn = np.abs(h_SR)**2

        h_DR_sq_mgn = np.abs(h_DR)**2 

        h_SD_sq_mgn = np.abs(h_SD)**2

        h_DD_tilde_sq_mgn = np.abs(h_DD_tilde)

        h_DW_sq_mgn = np.abs(h_DW)**2

        denominator = (2**r_P_bar - 1)

        first_portion = (1/h_DR_sq_mgn) * ((h_SR_sq_mgn * P_S/denominator) - No_R)

        second_portion = (1/h_DD_tilde_sq_mgn) * ((h_SD_sq_mgn * P_S / denominator) - No_D)

        third_portion = (No_W_unc **( 1 - 4 * err_min) - 1/No_W_unc) * (No_W / h_DW_sq_mgn )

        P_D_ub = np.array([first_portion,second_portion,third_portion,P_D_bar])
    
    P_D_opt = max(0,np.amin(P_D_ub))
    
    solutions = {}

    solutions["P_D"] = np.copy(P_D_opt)

    

    return solutions

def myf_algorithm_2(sys_param,channel,P_S_proportion_fixed):

    P_S = sys_param["P_S"]

    P_D_bar = sys_param["P_D_bar"]

    r_P_bar = sys_param["r_P_bar"]

    err_min = sys_param["err_min"]

    P_D_opt = min(P_S_proportion_fixed * P_S,P_D_bar)

    solutions = {}

    solutions["P_D"] = np.copy(P_D_opt)

    if (\
        (r_P_bar <= myf_r_P_R(sys_param, channel, solutions)) and \
       (r_P_bar <= myf_r_P_D(sys_param, channel, solutions)) and \
       (myf_DEP(sys_param, channel, solutions) >= err_min)\
       ):
       

       None
       
       
    else:
        
        P_D_opt = 0

    solutions["P_D"] = np.copy(P_D_opt)

    return solutions

         


