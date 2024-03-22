import numpy as np
# .....
# fParam.py
# .......
def myf_sys_param():

    sys_param = {}

    sys_param["BW"] = 20 * (10**6)

    sys_param["d_OS_range"] = 100

    sys_param["d_OD_range"] = 100

    sys_param["d_OR_range"] = 100

    sys_param["d_OW_range"] = 100

    sys_param["P_S_dBm"] = 23

    sys_param["P_S"] = 10 **(sys_param["P_S_dBm"]/10) / (10**3)

    sys_param["P_D_bar_dBm"] = 23

    sys_param["P_D_bar"] = 10 **(sys_param["P_D_bar_dBm"]/10) /(10 **3)

    sys_param["r_P_bar"] = 0.1

    sys_param["No_W_dBm"] = -160 #[dBm/Hz]

    sys_param["No_W"] = 10**(sys_param["No_W_dBm"] / 10)/(10**3) * sys_param["BW"] #[Watt]

    sys_param["No_W_unc_dB"] = 5

    sys_param["No_W_unc"] = 10**(sys_param["No_W_unc_dB"] / 10)

    sys_param["No_D_dBm"] = -160

    sys_param["No_D"] = 10**(sys_param["No_D_dBm"] / 10)/(10**3) * sys_param["BW"]

    sys_param["No_R_dBm"] = -160

    sys_param["No_R"] = 10**(sys_param["No_R_dBm"] / 10)/(10**3) * sys_param["BW"]

    sys_param["res_SI_dB"] = -100

    sys_param["res_SI"] = 10**(sys_param["res_SI_dB"] / 10)

    sys_param["err_min"] = 0.45

    sys_param["PL_exp"] = 3.5

    return sys_param



