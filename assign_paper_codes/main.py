import numpy as np

from fParam import *

from fChannel import *

from fAlgorithms import *

from fCalculations import *

from fPlot import *


sys_param = myf_sys_param()

Num_samples = 10000

x_axis_name = "No_W_unc_dB_cand" 

# x_axis_name = "P_S_dBm_cand"

if(x_axis_name == "P_S_dBm_cand"):
    x_axis_cand = np.arange(start=0 , stop=55, step=5)

elif(x_axis_name == "P_D_bar_dBm_cand"):

    x_axis_cand = np.arange(start=0 , stop=55, step=5)

elif(x_axis_name == "r_P_bar_cand"):

    x_axis_cand = np.arange(start=0 , stop=2.2, step=0.2)

elif(x_axis_name == "res_SI_dB_cand"):

    x_axis_cand = np.arange(start=-150 , stop=-40, step=10)

elif(x_axis_name == "err_min_cand"):

    x_axis_cand = np.arange(start=0 , stop=0.55, step=0.5)

elif(x_axis_name == "No_W_unc_dB_cand"):

    x_axis_cand = np.arange(start=0.5 , stop=5.5, step=0.5)

Num_algorithms = 5

r_C_R = np.zeros((Num_algorithms,np.size(x_axis_cand, axis =0 )))

r_P_R = np.zeros((Num_algorithms,np.size(x_axis_cand, axis =0 )))

r_P_D = np.zeros((Num_algorithms,np.size(x_axis_cand, axis =0 )))

DEP = np.zeros((Num_algorithms,np.size(x_axis_cand, axis =0 )))

Solutions_P_D = np.zeros((Num_algorithms,np.size(x_axis_cand, axis =0 )))

progress_per = 0

for ind1 in range(0,Num_samples):

    param_locations = {}

    param_locations["seed_seq"] = ind1

    locations = myf_locations(sys_param, param_locations)

    

    #Plot 

    # if(ind1==0): #Only at the first time

    #     myf_plot_locations(sys_param,locations)
     
    param_channel = {}

    param_channel["seed_seq"] = Num_samples + ind1
    
    channel = myf_channel(sys_param,locations,param_channel)


    r_C_R_temp = np.zeros((Num_algorithms, np.size(x_axis_cand,axis=0)))

    r_P_R_temp = np.zeros((Num_algorithms, np.size(x_axis_cand,axis=0)))

    r_P_D_temp = np.zeros((Num_algorithms, np.size(x_axis_cand,axis=0)))

    DEP_temp = np.zeros((Num_algorithms, np.size(x_axis_cand,axis=0)))

    Solutions_P_D_temp = np.zeros((Num_algorithms, np.size(x_axis_cand,axis=0)))

    for ind2 in range(0,np.size(x_axis_cand,axis=0)):

        if (x_axis_name == "P_S_dBm_cand" ):

            sys_param["P_S_dBm"] = x_axis_cand[ind2]

            sys_param["P_S"] = 10 **(sys_param["P_S_dBm"]/10)/(10**3)

               
        elif (x_axis_name =="P_D_bar_cand"):

            sys_param["P_D_bar_dBm"] = x_axis_cand[ind2]

            sys_param["P_D_bar"] = 10 **(sys_param["P_D_bar_dBm"]/10)/(10**3)

        elif (x_axis_name == "r_P_bar_cand"):

            sys_param["r_P_bar"] = x_axis_cand[ind2]

        elif (x_axis_name == "res_SI_dB_cand" ):

            sys_param["res_SI_dB"] = x_axis_cand[ind2]

            sys_param["res_SI"] = 10 **(sys_param["res_SI_dB"]/10)

            channel = myf_channel(sys_param, locations, param_channel)

        elif (x_axis_name =="err_min_cand"):

            sys_param["err_min"] = x_axis_cand[ind2]
        
        elif (x_axis_name == "No_W_unc_dB_cand"):
            
            sys_param["No_W_unc_dB"] = x_axis_cand[ind2]

            sys_param["No_W_unc"] = 10 **(sys_param["No_W_unc_dB"]/10)

        solutions_algorithm_1 = myf_algorihtm_1(sys_param,channel)  

        r_C_R_temp[0,ind2] = myf_r_C_R(sys_param,channel,solutions_algorithm_1)

        r_P_R_temp[0,ind2] = myf_r_P_R(sys_param,channel,solutions_algorithm_1)

        r_P_D_temp[0,ind2] = myf_r_P_D(sys_param,channel,solutions_algorithm_1)

        DEP_temp[0,ind2] = myf_DEP(sys_param,channel,solutions_algorithm_1)

        Solutions_P_D_temp[0,ind2]=solutions_algorithm_1["P_D"]

        

#         solutions_algorithm_2 = myf_algorithm_2(sys_param,channel,0.05)

        

#         r_C_R_temp[1,ind2] = myf_r_C_R(sys_param,channel,solutions_algorithm_2)

#         r_P_R_temp[1,ind2] = myf_r_P_R(sys_param,channel,solutions_algorithm_2)

#         r_P_D_temp[1,ind2] = myf_r_P_D(sys_param,channel,solutions_algorithm_2)

#         DEP_temp[1,ind2] = myf_DEP(sys_param,channel,solutions_algorithm_2)

#         Solutions_P_D_temp[1,ind2] = solutions_algorithm_2["P_D"]

#         solutions_algorithm_3 = myf_algorithm_2(sys_param,channel,0.01)

#         r_C_R_temp[2,ind2] = myf_r_C_R(sys_param,channel,solutions_algorithm_3)

#         r_P_R_temp[2,ind2] = myf_r_P_R(sys_param,channel,solutions_algorithm_3)

#         r_P_D_temp[2,ind2] = myf_r_P_D(sys_param,channel,solutions_algorithm_3)

#         Solutions_P_D_temp[2,ind2] = solutions_algorithm_3["P_D"]

#         solutions_algorithm_4 = myf_algorithm_2(sys_param,channel,0.001)

#         r_C_R_temp[3,ind2] = myf_r_C_R(sys_param,channel,solutions_algorithm_4)

#         r_P_R_temp[3,ind2] = myf_r_P_R(sys_param,channel,solutions_algorithm_4)

#         r_P_D_temp[3,ind2] = myf_r_P_D(sys_param,channel,solutions_algorithm_4)

#         DEP_temp[3,ind2] = myf_DEP(sys_param,channel,solutions_algorithm_4)

#         Solutions_P_D_temp[3,ind2] = solutions_algorithm_4["P_D"]

#         solutions_algorithm_5 = myf_algorithm_2(sys_param,channel,np.random.rand())

#         r_C_R_temp[4,ind2] = myf_r_C_R(sys_param,channel,solutions_algorithm_5)

#         r_P_R_temp[4,ind2] = myf_r_P_R(sys_param,channel,solutions_algorithm_5)

#         r_P_D_temp[4,ind2] = myf_r_P_D(sys_param,channel,solutions_algorithm_5)

#         Solutions_P_D_temp[4,ind2] = solutions_algorithm_5["P_D"]

#         progress = (ind1 * np.size(x_axis_cand,axis = 0) + ind2)\
#                     /(Num_samples * np.size(x_axis_cand,axis=0))*100
#         if(np.floor(progress/10) > progress_per):

#             progress_per = progress_per + 1

#             print("progress: ", progress)

#             print("....%d %% Done"%(progress_per*10))

#         ##average
        
        r_C_R = (ind1 + 1 - 1) / (ind1 + 1)*r_C_R + 1/(ind1 + 1) * r_C_R_temp

#         r_P_R = (ind1 + 1 - 1) / (ind1 + 1)*r_P_R + 1/(ind1 + 1) * r_P_R_temp

#         r_P_D = (ind1 + 1 - 1) / (ind1 + 1)*r_P_D + 1/(ind1 + 1) * r_P_D_temp

#         DEP = (ind1 + 1 - 1) / (ind1 + 1)*DEP + 1/(ind1 + 1) * DEP_temp

#         Solutions_P_D = (ind1 + 1 - 1) / (ind1 + 1)*Solutions_P_D + 1/(ind1 + 1) * Solutions_P_D_temp

print("r_C_R(alogrithm 1):", r_C_R[0,:])

# print("Solutions_P_D(algorihtm 1): ", Solutions_P_D[0,:])

# print("r_C_R(alogrithm 2):", r_C_R[1,:])

# print("r_C_R(alogrithm 3):", r_C_R[2,:])

# print("r_C_R(alogrithm 4):", r_C_R[3,:])

# print("r_C_R(alogrithm 5):", r_C_R[4,:])

# myf_plot_r_C_R(sys_param,x_axis_cand,r_C_R,r_P_R,r_P_D,x_axis_name)

# myf_plot_Solutions_P_D(sys_param,x_axis_cand,Solutions_P_D,x_axis_name)

# myf_plot_DEP(sys_param,x_axis_cand,DEP,x_axis_name)

# print("Finished")