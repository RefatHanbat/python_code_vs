import time
import numpy as np

#fCalculations.py

def myf_r_C_R(sys_param, channel, solutions):

    No_R = sys_param["No_R"]

    h_DR = channel["h_DR"]

    P_D = solutions["P_D"]

    numerator = np.abs(h_DR)**2 * P_D

    denominator = No_R

    r_C_R = np.log2(1 + numerator/denominator)

    return r_C_R
def myf_r_P_R(sys_param, channel, solutions):

    P_S = sys_param["P_S"]

    No_R = sys_param["No_R"]

    h_SR = sys_param["h_SR"]

    h_DR = sys_param["h_DR"]

    P_D = solutions["P_D"]

    numerator = np.abs(h_SR)**2*P_S

    denominator = np.abs(h_DR)**2 * P_D + No_R

    r_P_R = np.log2(1 + numerator/denominator)

    return r_P_R

def myf_r_P_D(sys_param, channel, solutions):

    P_S = sys_param["P_S"]

    No_D = sys_param["No_D"]

    h_SD = channel["h_SD"]

    h_DD_tilde = channel["h_DD_tilde"]

    P_D = solutions["P_D"]

    numerator = np.abs(h_SD)**2*P_S

    denominator = np.abs(h_DD_tilde)**2*P_D + No_D

    r_P_D = np.log2(1 + numerator/denominator)

    return r_P_D

def myf_DEP(sys_param,channel,solutions):

    No_W = sys_param["No_W"]

    No_W_unc = sys_param["No_W_unce"]

    h_DW = channel["h_DW"]

    P_D = solutions["P_D"]

    numerator = np.abs(h_DW)**2*P_D + (1/No_W_unc) * No_W

    denominator = (1/No_W_unc)*No_W

    DEP = (1/2)*(1-(1/2*np.log(No_W_unc))) * np.log(numerator/denominator)