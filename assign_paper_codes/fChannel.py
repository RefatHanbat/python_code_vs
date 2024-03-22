import numpy as np

# fChannel.py

def myf_locations(sys_param,param_locations):

    d_OS_range = sys_param["d_OS_range"]
    d_OD_range = sys_param["d_OD_range"]
    d_OR_range = sys_param["d_OR_range"]
    d_OW_range = sys_param["d_OW_range"]
    seed_seq = param_locations["seed_seq"]

    rng = np.random.default_rng(seed=seed_seq)

    locations = {}

    locations["d_OS"] = 1 + (d_OS_range - 1)

    locations["d_OD"] = 1 + (d_OD_range - 1)

    locations["d_OR"] = 1 + (d_OR_range - 1)

    locations["d_OW"] = 1 + (d_OW_range - 1)

    locations["angle_OS"] = np.pi

    locations["angle_OD"] = 0

    locations["angle_OR"] = (1/2)*np.pi

    locations["angle_OW"] = (3/2) * np.pi

    locations["positions_S"] = np.multiply(np.vstack((locations["d_OS"], locations["d_OS"])), \
                                       np.vstack((np.cos(locations["angle_OS"]), \
                                                  np.sin(locations["angle_OS"]))))
    
    locations["positions_D"] = np.multiply(np.vstack((locations["d_OD"], locations["d_OD"])), \
                                       np.vstack((np.cos(locations["angle_OD"]), \
                                                  np.sin(locations["angle_OD"]))))
    
    locations["positions_R"] = np.multiply(np.vstack((locations["d_OR"], locations["d_OR"])), \
                                       np.vstack((np.cos(locations["angle_OR"]), \
                                                  np.sin(locations["angle_OR"]))))
    
    locations["positions_W"] = np.multiply(np.vstack((locations["d_OW"], locations["d_OW"])), \
                                       np.vstack((np.cos(locations["angle_OW"]), \
                                                  np.sin(locations["angle_OW"]))))


    return locations

def myf_channel(sys_param, locations, param_channel):

    res_SI = sys_param["res_SI"]

    PL_exp = sys_param["PL_exp"]

    positions_S = locations["positions_S"]

    positions_D = locations["positions_D"]

    positions_R = locations["positions_R"]

    positions_W = locations["positions_W"]

    seed_seq = param_channel["seed_seq"]


    rng = np.random.default_rng(seed=seed_seq)

    channel = {}

    # PL_SD = 10**(-3) * np.sqrt(np.sum(np.power(positions_D[:,0]-positions_S[:,0], 2))**(-PL_exp))

    PL_SD = 10**(-3) * np.sqrt(np.sum(np.power(positions_D[:,0]-positions_S[:,0],2),0))**(-PL_exp)

    # PL_SR = 10**(-3) * np.sqrt(np.sum(np.power(positions_R[:,0]-positions_S[:,0], 2)))**(-PL_exp)

    PL_SR = 10**(-3) * np.sqrt(np.sum(np.power(positions_R[:,0]-positions_S[:,0],2),0))**(-PL_exp)

    PL_SW = 10**(-3) * np.sqrt(np.sum(np.power(positions_W[:,0]-positions_S[:,0],2),0))**(-PL_exp)

    # PL_SW = 10**(-3) * np.sqrt(np.sum(np.power(positions_W[:,0]-positions_S[:,0],2)))**(-PL_exp)

    PL_DR = 10**(-3) * np.sqrt(np.sum(np.power(positions_R[:,0]-positions_D[:,0],2),0))**(-PL_exp)

    # PL_DR = 10**(-3) * np.sqrt(np.sum(np.power(positions_R[:,0]-positions_D[:,0], 2)))**(-PL_exp)

    PL_DW = 10**(-3) * np.sqrt(np.sum(np.power(positions_W[:,0]-positions_D[:,0],2),0))**(-PL_exp)

    # PL_DW = 10**(-3) * np.sqrt(np.sum(np.power(positions_W[:,0]-positions_D[:,0], 2)))**(-PL_exp)


    PL_DD = res_SI

    h_SD = rng.normal(loc=0, scale=np.sqrt(1/2)) + 1j * rng.normal(loc=0, scale=np.sqrt(1/2))

    h_SR = rng.normal(loc=0, scale=np.sqrt(1/2)) + 1j * rng.normal(loc=0, scale=np.sqrt(1/2))

    h_SW = rng.normal(loc=0, scale=np.sqrt(1/2)) + 1j * rng.normal(loc=0, scale=np.sqrt(1/2))

    h_DR = rng.normal(loc=0, scale=np.sqrt(1/2)) + 1j * rng.normal(loc=0, scale=np.sqrt(1/2))

    h_DW = rng.normal(loc=0, scale=np.sqrt(1/2)) + 1j * rng.normal(loc=0, scale=np.sqrt(1/2))

    h_DD_tilde = rng.normal(loc=0, scale=np.sqrt(1/2)) + 1j * rng.normal(loc=0, scale=np.sqrt(1/2))

    channel["h_SD"] = np.sqrt(PL_SD) * h_SD

    channel["h_SR"] = np.sqrt(PL_SR) * h_SR

    channel["h_SW"] = np.sqrt(PL_SW) * h_SW

    channel["h_DR"] = np.sqrt(PL_DR) * h_DR

    channel["h_DW"] = np.sqrt(PL_DW) * h_DW

    channel["h_DD_tilde"] = np.sqrt(PL_DD) * h_DD_tilde

    return channel


