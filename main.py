import numpy as np

import matplotlib.pyplot as plt

N = 1000

ip = np.random.rand(N) > 0.5

ip = ip.astype(int)

s = 2 * ip - 1

nTx = 2

nRx = 2

Es = 1

eb_no_db = np.arange(1)

def design_precoder(Ns,optimal_Q,V,check_csit):

    type = check_csit

    if(type=="without_csit"):

        P = np.sqrt(Es / Ns) * np.vstack([
               np.eye(Ns),                
             np.zeros((nTx-Ns, Ns))           
          ])
    elif(type=="svd"):

        Ns = len(optimal_Q)

        q_values = np.array(list(optimal_Q.values()))

        q_sqrt = np.sqrt(q_values)

        q_sqrt = np.pad(q_sqrt, (0, nTx - len(q_sqrt)))

        Q_matrix = np.diag(q_sqrt)

        P = V[:,0:Ns] @ Q_matrix[0:Ns,0:Ns]

    return P





def findOptimalQ(noise_variance,S):

    max_length_q_values = {}

    for ii in range(1, len(S)+1):
        
        subset = S[:ii]

        summation_part = 0

        q_values = {}

        break_outer = False

        for kk in range(len(subset)):

            summation_part = summation_part + (noise_variance/np.square(subset[kk]))

            mu = (1/len(subset)) * (Es + summation_part)
        
        for mm in range(len(subset)):

            q_opt = mu - (noise_variance/np.square(subset[mm]))

            if (q_opt > 0):

                q_values[mm] = q_opt

            else:

                break_outer = True

                break

        if break_outer:

            break
        if len(q_values) > len(max_length_q_values):

            max_length_q_values = q_values

    return max_length_q_values
    
        

def water_filling():

    for jj in range(len(eb_no_db)):

        snr = 10**(eb_no_db[jj]/10)

        noise_variance = 1/(snr)

        Rz = noise_variance * np.eye(nRx)
        
        Rz_inverse = np.linalg.inv(Rz)

        modulated_effective_bits = []

        demodulated_effective_bits = []

        equalize_bits = []



        for ii in range(N//2):

            H = np.sqrt(1/(2)) * (np.random.randn(nRx, nTx) + 1j * np.random.randn(nRx, nTx))

            H_hermitian = np.conjugate(H)

            H_hermitian = np.transpose(H_hermitian)

            Ns = np.linalg.matrix_rank(H)

            U,S,Vh = np.linalg.svd(H,full_matrices=False)

            V = np.conjugate(Vh)

            V = np.transpose(V)

            optimal_Q = findOptimalQ(noise_variance,S)

            # precoder = design_precoder(Ns,optimal_Q,V,"without_csit")

            P = design_precoder(Ns,optimal_Q,V,"svd")

            P_hermitian = np.conjugate(P)

            P_hermitian = np.transpose(P_hermitian)

            n = (1/np.sqrt(2)) * (np.random.randn(2, 1) + 1j * np.random.randn(2, 1))

            x1 = s[2 * ii]
                              
            x2 = s[2 * ii + 1]

            x = np.array([[x1], [x2]])

            if(P.shape == (2,2)):

                modulated_effective_bits.append(x1)

                modulated_effective_bits.append(x2)

                y = H @ P @ x + n * (10 ** (-eb_no_db[jj] / 20))

                W = P_hermitian @ H_hermitian @ Rz_inverse

                x_hat = W @ y

                x_hat = np.real(x_hat > 0).astype(int)

                equalize_bits.append(x_hat[0, 0])

                equalize_bits.append(x_hat[1, 0])
                

            elif(P.shape == (2,1)):

                modulated_effective_bits.append(x1)

                x = np.array([[x1]])

                y = H @ P @ x + n * (10 ** (-eb_no_db[jj] / 20))

                W = P_hermitian @ H_hermitian @ Rz_inverse

                x_hat = W @ y

                x_hat = np.real(x_hat > 0).astype(int)
            
                equalize_bits.append(x_hat[0, 0])
        
        for ff in range(len(modulated_effective_bits)):

            print(modulated_effective_bits[ff])

            if(modulated_effective_bits[ff] < 0 ):

                demodulated_effective_bits.append(0)

            else:

                demodulated_effective_bits.append(1)

        print(demodulated_effective_bits)

        print(equalize_bits)

        nErr = np.sum(demodulated_effective_bits != equalize_bits)

        print(nErr)
        

        


            

            
water_filling()            


