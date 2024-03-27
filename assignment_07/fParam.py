import numpy as np

def myf_sys_param():

    sys_param = {}

    sys_param["Nt"] = 2

    sys_param["Nr"] = 2

    sys_param["Ns"] = min(sys_param["Nt"],sys_param["Nr"])

    sys_param["Es"] = 1

    sys_param["EsoverNO_dB"] = 5

    sys_param["EsoverNO"] = 10**(sys_param["EsoverNO_dB"]/10)

    sys_param["NO"] = sys_param["Es"]/(sys_param["EsoverNO"])

    sys_param["constellation_type"] = "QPSK"

    if(sys_param["constellation_type"] == "BPSK"):

        sys_param["constellation_size"] =2

    elif(sys_param["constellation_type"] == "QPSK"):

        sys_param["constellation_size"] = 4

    return sys_param