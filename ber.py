import numpy as np

import matplotlib.pyplot as plt

N = 20000

seed_seq = np.random.randint(0, N)

rng = np.random.default_rng(seed=seed_seq)

ip = np.random.rand(N) > 0.5

ip = ip.astype(int)

s = 2 * ip - 1

nTx = 2

nRx = 2

Es = 1

eb_no_db = np.arange(0,25)

def design_precoder(Ns,optimal_Q,V,check_csit):

    type = check_csit

    if(type=="without_csit"):

        P = np.sqrt(Es / Ns) * np.vstack([
               np.eye(Ns),                
             np.zeros((nTx-Ns, Ns))           
          ])
    elif(type=="svd"):

        Ns = np.count_nonzero(np.diag(optimal_Q))

        Q_matrix = np.sqrt(optimal_Q)

        

        P = V[:,0:Ns] @ Q_matrix[0:Ns,0:Ns]

    return P,Ns





def findOptimalQ_origin(noise_variance,S):

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


def findOptimalQ(No_var, S):

    Ns = min(nTx, nRx)

    S_mtr = np.diag(S)

    Q_mtr_opt = np.zeros((Ns, Ns), dtype=np.float32)

    C_opt = 0

    # loop for Ns number of stream
    for ind1 in range(0, Ns):
        flag_some_q_neg = 0
        Q_mtr_temp = np.zeros((Ns, Ns), dtype=np.float32)
        Ns_temp = ind1 + 1
        C_temp = 0
        sum_temp = 0

        # loop for ind2 is double than when ind1
        for ind2 in range(0, Ns_temp):
            sum_temp = sum_temp + No_var / S_mtr[ind2][ind2]

        # Getting Mu values
        mu_temp = (1 / Ns_temp) * (Es + sum_temp)

        for ind2 in range(0, Ns_temp):
            Q_mtr_temp[ind2][ind2] = mu_temp - No_var / S_mtr[ind2][ind2]

            if Q_mtr_temp[ind2][ind2] < 0:
                flag_some_q_neg = 1
                break

        if flag_some_q_neg != 1:
            for ind2 in range(0, Ns_temp):
                C_temp = C_temp + np.log2(1 + (1 / No_var) * (S_mtr[ind2][ind2] ** 2) * Q_mtr_temp[ind2][ind2])

            if C_temp > C_opt:
                Q_mtr_opt = np.copy(Q_mtr_temp)
                Ns_opt = np.copy(Ns_temp)
                C_opt = np.copy(C_temp)
                
    return Q_mtr_opt

def Rx_svd(type):

    print(f"Progress: Completed simulation for {type} equalization")
    equalization_name = type

    nErr = []
    
    for jj in range(len(eb_no_db)):

        snr = 10**(eb_no_db[jj]/10)

        noise_variance = 1/(snr)

        Rz = noise_variance * np.eye(nRx)

        Rz_inverse = np.linalg.inv(Rz)

        ii = 0

        equalize_bits = []

        while ii < len(s):

            H = rng.normal(loc=0, scale=np.sqrt(1/2), size=(nTx, nRx)) + 1j * rng.normal(loc=0, scale=np.sqrt(1/2), size=(nTx, nRx))

            H_hermitian = np.conjugate(H)

            H_hermitian = np.transpose(H_hermitian)

            U, S, Vh = np.linalg.svd(H, full_matrices=True)

            V = np.conjugate(Vh)

            V = np.transpose(V)

            n = (1/np.sqrt(2)) * (np.random.randn(2, 1) + 1j * np.random.randn(2, 1))

            rank = np.linalg.matrix_rank(H)

            optimal_Q = findOptimalQ(noise_variance, S)

            P, Ns = design_precoder(rank, optimal_Q, V, "svd")

            P_hermitian = np.conjugate(P)

            P_hermitian = np.transpose(P_hermitian)

            if P.shape == (2, 1):

                x1 = s[ii]

                x = np.array([[x1]])

                y = H @ P @ x + n * (10 ** (-eb_no_db[jj] / 20))

                if equalization_name == "MF":

                    W = P_hermitian @ H_hermitian @ Rz_inverse

                elif equalization_name == "ZF":

                    W = np.linalg.inv(P_hermitian @ H_hermitian @ Rz_inverse @ H @ P) @ P_hermitian @ H_hermitian @ Rz_inverse

                elif equalization_name == "MMSE":

                    W = np.linalg.inv(np.eye(Ns) + P_hermitian @ H_hermitian @ Rz_inverse @ H @ P) @ P_hermitian @ H_hermitian @ Rz_inverse

                x_hat = np.dot(W, y)

                x_hat = np.real(x_hat > 0).astype(int)

                equalize_bits.append(x_hat[0, 0])

                ii = ii + 1

            elif P.shape == (2, 2):

                if (len(s) - ii) < 2:

                    x1 = s[ii]

                    x = np.array([[x1], [0]])

                    y = H @ P @ x + n * (10 ** (-eb_no_db[jj] / 20))

                    if equalization_name == "MF":

                        W = P_hermitian @ H_hermitian @ Rz_inverse

                    elif equalization_name == "ZF":

                        W = np.linalg.inv(P_hermitian @ H_hermitian @ Rz_inverse @ H @ P) @ P_hermitian @ H_hermitian @ Rz_inverse
                        
                    elif equalization_name == "MMSE":

                        W = np.linalg.inv(np.eye(Ns) + P_hermitian @ H_hermitian @ Rz_inverse @ H @ P) @ P_hermitian @ H_hermitian @ Rz_inverse

                    x_hat = np.dot(W, y)

                    x_hat = np.real(x_hat > 0).astype(int)

                    equalize_bits.append(x_hat[0, 0])

                    ii = ii + 2

                else:
                    x1 = s[ii]

                    x2 = s[ii + 1]

                    x = np.array([[x1], [x2]])

                    y = H @ P @ x + n * (10 ** (-eb_no_db[jj] / 20))

                    if equalization_name == "MF":

                        W = P_hermitian @ H_hermitian @ Rz_inverse

                    elif equalization_name == "ZF":

                        W = np.linalg.inv(P_hermitian @ H_hermitian @ Rz_inverse @ H @ P) @ P_hermitian @ H_hermitian @ Rz_inverse

                    elif equalization_name == "MMSE":

                        W = np.linalg.inv(np.eye(Ns) + P_hermitian @ H_hermitian @ Rz_inverse @ H @ P) @ P_hermitian @ H_hermitian @ Rz_inverse

                    x_hat = np.dot(W, y)

                    x_hat = np.real(x_hat > 0).astype(int)

                    equalize_bits.append(x_hat[0, 0])

                    equalize_bits.append(x_hat[1, 0])

                    ii = ii + 2

        nErr.append(np.sum(ip != equalize_bits))

    return nErr

        

mf_svd = Rx_svd("MF")

zf_svd = Rx_svd("ZF")

mmse_svd = Rx_svd("MMSE")

plt.figure()

plt.semilogy(eb_no_db, mf_svd, 'c^-', linewidth=2, markersize=8, label='nTx=2 nRx=2 (mimo w/o csit RxMF)(EPA SVD P)')

plt.semilogy(eb_no_db, zf_svd, 'k>-', linewidth=2, markersize=8, label='nTx=2 nRx=2 (mimo w/o csit RxZF)(EPA SVD P)')

plt.semilogy(eb_no_db, mmse_svd, 'y<-', linewidth=2, markersize=8, label='nTx=2 nRx=2 (mimo w/o csit RxMMSE)(EPA SVD P)')

plt.grid(True)
plt.legend()
plt.xlabel('Eb/No, dB')
plt.ylabel('Bit Error Rate')
plt.title('SNR-dB')
plt.show()






