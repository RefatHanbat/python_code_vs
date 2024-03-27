import numpy as np

import copy

from fParam import*

from fChannel import*

from fAlgorithms import*

from fCalculations import*

from fPlot import*

sys_param = myf_sys_param()


Num_samples = 1000

x_axis_name = "EsoverNO_dB_cand"

if(x_axis_name == "EsoverNO_dB_cand"):

    x_axis_cand = np.arange(start=0, stop=30, step=5)


elif(x_axis_name == "EboverNO_dB_cand"):

    x_axis_cand = np.arange(start=0, stop=30, step=5)

Num_algorihtms = 7

###constellation symbols and bits (Unit average symbol energy)

constellation = myf_constellation(sys_param)

Num_errors_symbol = np.zeros((Num_algorihtms, np.size(x_axis_cand, axis=0)))

Num_errors_bit = np.zeros((Num_algorihtms, np.size(x_axis_cand, axis=0)))

Num_symbols_eff = np.zeros((Num_algorihtms, np.size(x_axis_cand, axis=0)))

Num_bits_eff = np.zeros((Num_algorihtms, np.size(x_axis_cand, axis=0)))

progress_per = 0

for ind1 in range(0,Num_samples):

    ##channel 
    
    param_channel = {}

    param_channel["seed_seq"] = Num_samples + ind1

    channel = myf_channel(sys_param, param_channel)

    

    ##Step1 : Generate Ns random indices of the constellation

    param_symbol_indices = {}

    param_symbol_indices["seed_seq"] = Num_samples + 1 + ind1


    symbol_indices = myf_symbol_indices(sys_param, param_symbol_indices)


    x_vec = myf_x_vec(sys_param,constellation, symbol_indices)

    ##step(2):Encoding i.e placing the constellation symbols ontx x_vec according to 
    ##symbol_indices(note R_x_vec= I)

    Num_errors_symbol_temp = np.zeros((Num_algorihtms, np.size(x_axis_cand, axis=0)))

    Num_errors_bit_temp = np.zeros((Num_algorihtms, np.size(x_axis_cand, axis=0)))

    Num_symbols_eff_temp = np.zeros((Num_algorihtms, np.size(x_axis_cand, axis=0)))

    Num_bits_eff_temp = np.zeros((Num_algorihtms, np.size(x_axis_cand, axis=0)))

    ##algorithms###
    
    for ind2 in range(0, np.size(x_axis_cand, axis=0)):

        if(x_axis_name == "EsoverNO_dB_cand"):

            sys_param["EsoverNO_dB"] = x_axis_cand[ind2]

            sys_param["EsoverNO"] = 10**(sys_param["EsoverNO_dB"]/10)

            sys_param["NO"] = sys_param["Es"] / sys_param["EsoverNO"]

        elif(x_axis_name=="EboverNO_dB_cand"):
            
            sys_param["EsoverNO_dB"]=x_axis_cand[ind2] \
                + 10*np.log10(np.log2(sys_param["constellation_size"]))
            

            sys_param["EsoverNO"] = 10**(sys_param["EsoverNO_dB"]/ 10)

            sys_param["NO"] = sys_param["Es"] / sys_param["EsoverNO"]
            
        ###Equal power allocation (EPA) precoder ###
        
        P_mat , P_mat_reduced, Ns_opt = myf_P(sys_param, channel, mode="EPA")


        ##step(3): Reduce symbols if Ns_opt < Ns

        x_vec_reduced = x_vec[0:Ns_opt,:]

#         ##print("x_vec_reduced:",x_vec_reduced)

#         ##step(4):Precoding i.e s_vec = P_mat@x_vec (note tre(R_s_vec) = Es)

        s_vec = P_mat_reduced@x_vec_reduced

        ##Step(5):Received symbol vector i.e , y_vec = H_mat * s_vec + z_vec

        param_y_vec = {}

        param_y_vec["seed_seq"] = Num_samples + 2 + ind1

        param_y_vec["s_vec"] = s_vec

        y_vec = myf_y_vec(sys_param, channel, param_y_vec)

        ##Algorithm 1: MIMO w/o CSIT - RxMF(EPA P)

        solutions_algorithm_1 = myf_algorithm_1(sys_param, channel, P_mat_reduced, Ns_opt, mode="RxMF")

        param_Num_errors = {}

        param_Num_errors["symbol_indices"] = symbol_indices

        param_Num_errors["Ns_opt"] = Ns_opt

        param_Num_errors["y_vec"] = y_vec

        Num_errors = myf_Num_errors(sys_param,constellation,solutions_algorithm_1,param_Num_errors)

        Num_errors_symbol_temp[0][ind2] = Num_errors["Num_errors_symbol"]

        Num_errors_bit_temp[0][ind2] = Num_errors["Num_errors_bit"]

        Num_symbols_eff_temp[0][ind2] = Ns_opt

        Num_bits_eff_temp[0][ind2] = Ns_opt * np.log2(sys_param["constellation_size"])

        ##Algorithm 2: MIMO w/o CSIT - RxZF (EPA P)

        solutions_algorithm_2 = myf_algorithm_1(sys_param, channel, P_mat_reduced, Ns_opt, mode="RxZF")

        param_Num_errors = {}

        param_Num_errors["symbol_indices"] = symbol_indices

        param_Num_errors["Ns_opt"] = Ns_opt

        param_Num_errors["y_vec"] = y_vec

        Num_errors = myf_Num_errors(sys_param,constellation,solutions_algorithm_2,param_Num_errors)

        Num_errors_symbol_temp[1][ind2] = Num_errors["Num_errors_symbol"]

        Num_errors_bit_temp[1][ind2] = Num_errors["Num_errors_bit"]

        Num_symbols_eff_temp[1][ind2] = Ns_opt

        Num_bits_eff_temp[1][ind2] = Ns_opt * np.log2(sys_param["constellation_size"])

        ##Algorithm 3: MIMO w/o CSIT -RxMMSE(EPA P)


        solutions_algorithm_3 = myf_algorithm_1(sys_param,channel, P_mat_reduced, Ns_opt,mode="RxMMSE")

        Num_errors = myf_Num_errors(sys_param,constellation, solutions_algorithm_3, param_Num_errors)

        Num_errors_symbol_temp[2][ind2] = Num_errors["Num_errors_symbol"]
        
        Num_errors_bit_temp[2][ind2] = Num_errors["Num_errors_bit"]

        Num_symbols_eff_temp[2][ind2] = Ns_opt

        Num_bits_eff_temp[2][ind2] = Ns_opt * np.log2(sys_param["constellation_size"])

        ##svd precoder##

        P_mat , P_mat_reduced, Ns_opt = myf_P(sys_param, channel, mode="SVD")

        x_vec_reduced = x_vec[0:Ns_opt,:]

        s_vec = P_mat_reduced@x_vec_reduced

        param_y_vec = {}

        param_y_vec["seed_seq"] = Num_samples + 2 + ind1 

        param_y_vec["s_vec"] = s_vec

        y_vec = myf_y_vec(sys_param,channel,param_y_vec)

        ##Algorithm 4: MIMO w/o CSIT - RxMF(SVD P)

        solutions_algorithm_4 = myf_algorithm_1(sys_param,channel,P_mat_reduced,Ns_opt,mode="RxMF")

        param_Num_errors = {}

        param_Num_errors["symbol_indices"] = symbol_indices

        param_Num_errors["Ns_opt"] = Ns_opt

        param_Num_errors["y_vec"] = y_vec

        Num_errors = myf_Num_errors(sys_param,constellation,solutions_algorithm_4,param_Num_errors)

        Num_errors_symbol_temp[3][ind2] = Num_errors["Num_errors_symbol"]

        Num_errors_bit_temp[3][ind2] = Num_errors["Num_errors_bit"]

        Num_symbols_eff_temp[3][ind2] = Ns_opt

        Num_bits_eff_temp[3][ind2] = Ns_opt * np.log2(sys_param["constellation_size"])

        ##Algorithm 5: MIMO w/o CSIT -RxZF (SVD P)

        solutions_algorithm_5 = myf_algorithm_1(sys_param,channel,P_mat_reduced,Ns_opt,mode="RxZF")

        param_Num_errors = {}

        param_Num_errors["symbol_indices"] = symbol_indices

        param_Num_errors["Ns_opt"] = Ns_opt

        param_Num_errors["y_vec"] = y_vec

        Num_errors = myf_Num_errors(sys_param,constellation,solutions_algorithm_5,param_Num_errors)

        Num_errors_symbol_temp[4][ind2] = Num_errors["Num_errors_symbol"]

        Num_errors_bit_temp[4][ind2] = Num_errors["Num_errors_bit"]

        Num_symbols_eff_temp[4][ind2] = Ns_opt

        Num_bits_eff_temp[4][ind2] = Ns_opt * np.log2(sys_param["constellation_size"])

        ##Algorithm 6: MIMO w/o CSIT -RxMMSE (SVD P)

        solutions_algorithm_6 = myf_algorithm_1(sys_param,channel,P_mat_reduced,Ns_opt,mode="RxZF")

        param_Num_errors = {}

        param_Num_errors["symbol_indices"] = symbol_indices

        param_Num_errors["Ns_opt"] = Ns_opt

        param_Num_errors["y_vec"] = y_vec

        Num_errors = myf_Num_errors(sys_param,constellation,solutions_algorithm_6,param_Num_errors)

        Num_errors_symbol_temp[5][ind2] = Num_errors["Num_errors_symbol"]

        Num_errors_bit_temp[5][ind2] = Num_errors["Num_errors_bit"]

        Num_symbols_eff_temp[5][ind2] = Ns_opt

        Num_bits_eff_temp[5][ind2] = Ns_opt * np.log2(sys_param["constellation_size"])

        ##Algorithm 7 : MIMO w/ CSIT (SVD P)

        solutions_algorithm_7 = myf_algorithm_2(sys_param, channel, mode="SVD")

        Num_errors = myf_Num_errors(sys_param,constellation,solutions_algorithm_7,param_Num_errors)

        Num_errors_symbol_temp[6][ind2] = Num_errors["Num_errors_symbol"]

        Num_errors_bit_temp[6][ind2] = Num_errors["Num_errors_bit"]

        Num_symbols_eff_temp[6][ind2] = Ns_opt

        Num_bits_eff_temp[6][ind2] = Ns_opt * np.log2(sys_param["constellation_size"])
    
        ##progress##

        progress = (ind1*np.size(x_axis_cand, axis=0) + ind2)\
                  /(Num_samples*np.size(x_axis_cand,axis=0))*100
        
        if(np.floor(progress/10) > progress_per):

            progress_per =  progress_per + 1

            print("Progress: ", progress)

            print("...%d %%Done..."%(progress_per*10))
    ##Accumulation
         
    Num_errors_symbol = Num_errors_symbol + Num_errors_symbol_temp

    Num_errors_bit = Num_errors_bit + Num_errors_bit_temp

    Num_symbols_eff = Num_symbols_eff + Num_symbols_eff_temp

    Num_bits_eff = Num_bits_eff + Num_bits_eff_temp

SER = np.divide(Num_errors_symbol,Num_symbols_eff)

BER = np.divide(Num_errors_bit,Num_bits_eff)


# print("Num_symbols_eff: ",Num_symbols_eff )

# print("Num_bits_eff : ",Num_bits_eff )

# print("BER (MIMO w/o CSIT - RxMF (EPA P)) : ", BER[0, :])

# print("BER (MIMO w/o CSIT - RxZF (EPA P)) : ", BER[1, :])

# print("BER (MIMO w/o CSIT - RxMMSE (EPA P)) : ", BER[2, :])

# print("BER (MIMO w/o CSIT - RxMF (SVD P)) : ", BER[3, :])

# print("BER (MIMO w/o CSIT - RxZF (SVD P)) : ", BER[4, :])

# print("BER (MIMO w/o CSIT - RxMMSE (SVD P)) : ", BER[5, :])

# print("BER (MIMO w/ CSIT -SVD P) : ", BER[5, :])

myf_plot_BER(sys_param,x_axis_cand,BER,x_axis_name)

print("Finished")





    

    
