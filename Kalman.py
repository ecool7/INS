import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np
R = 6371000
g = 9.81
k = 5000
T = 1
nu_k = 0.8
A = 0.001
beta = 1.13e-11
w_dr = 4.8*10e-8
#####      МАТРИЦЫ  #########
dV = np.zeros(k)
Fn = np.zeros(k)
wdr = np.zeros(k)
Zkdv1 = np.zeros(k)
Zkdv2 = np.zeros(k)
Zkdv3 = np.zeros(k)
Zkdv4 = np.zeros(k)
Zkdv5 = np.zeros(k)
Zkdv6 = np.zeros(k)
Zkdv7 = np.zeros(k)
Zkdv8 = np.zeros(k)
Zkdv9 = np.zeros(k)
F_MATR = np.array([[1,-g*T,0],[T/R, 1,T],[0,0,1-beta*T]])

for i in range(1,k):

    b = np.array([dV[i-1],Fn[i-1],wdr[i-1]])

    a = np.array([0, 0, A*((2*beta)**(0.5))*(random.uniform(-1,1))])+ np.dot(F_MATR, b)
    dV[i] = a[0]
    Fn[i] = a[1]
    wdr[i] = a[2]
    Zkdv1[i] = dV[i] - (np.sqrt(1))*(random.uniform(-1,1))
    Zkdv2[i]= dV[i] - (np.sqrt(0.01))*(random.uniform(-1,1))
    Zkdv3[i] = dV[i] - (np.sqrt(0.0001))*(random.uniform(-1,1))



np.savetxt("file_dV.txt", dV)
np.savetxt("file_Fn.txt", Fn)
np.savetxt("file_wdr.txt", wdr)
np.savetxt("file_Zkdv1.txt", Zkdv1)
np.savetxt("file_Zkdv2.txt", Zkdv2)
np.savetxt("file_Zkdv3.txt", Zkdv3)




Z = np.loadtxt("file_Zkdv1.txt")
dV = np.loadtxt("file_dV.txt")
Fn = np.loadtxt("file_Fn.txt")
wdr = np.loadtxt("file_wdr.txt")
Z1= np.loadtxt("file_Zkdv1.txt")
Z2 = np.loadtxt("file_Zkdv2.txt")
Z3 = np.loadtxt("file_Zkdv3.txt")



# np.loadtxt("file_Zkdv2.txt", Zkdv2)
# np.loadtxt("file_Zkdv3.txt", Zkdv3)

R = 1
Q = 1E-16
X0 = np.zeros((k,3))
H = np.array([1, 0, 0])
G = np.array([[0,0,0],[0, 0,0],[0,0,1]])
Xk = X0
I = np.eye(3)
A = np.array([[0,g,0],[-1/R, 0,1],[0,0,-beta]])
F = I+A*T+A**(2)*T**(2/2)
est_X0 = np.array([0, 0, 0])
est_X_prev = est_X0
P0 = np.array([[1,0,0],[0, 1,0],[0,0,1]])
P_prev = P0
Nuk = np.array([0, 0, 0])
#рикати
prcnt10 = k/100
p_count = 0

X_est = np.zeros((k,3)) # с домиком - оценка
Pk = np.zeros((3,3)) # 3x3
Pk_k_1 = np.zeros((3,3)) # 3x3 Pk/k-1

T=np.arange(1,k)
NN_N=200
PP_1 = np.zeros(NN_N)
PP_2 = np.zeros(NN_N)
PP_3 = np.zeros(NN_N)



Q_fsk = 1e-13
X_est_fsk = np.zeros((k, 3))  # с домиком - оценка
Pk_fsk = np.diag([10, 1e-6, 1e-9])  # 3x3
Pk_k_1_fsk = np.zeros((3, 3))  # 3x3 Pk/k-1
KK = np.dot(Pk_k_1_fsk, np.transpose(H)) / np.transpose(np.dot(H, np.dot(Pk_k_1_fsk, np.transpose(H))) + R)  # усиление
for i in range(1, k):
    Pk_k_1_fsk = np.dot(F, np.dot(Pk_fsk, np.transpose(F))) + np.dot(G, Q * G)
    KK = np.dot(Pk_k_1_fsk, H[:, np.newaxis]) / (np.dot(H, np.dot(Pk_k_1_fsk, H[:, np.newaxis])) + R)
    X_old_fsk = X_est_fsk[i - 1]
    A = np.dot(H, np.dot(F, X_old_fsk[:, np.newaxis]))
    KKK = np.array([KK[0, 0], KK[1, 0], KK[2, 0]])
    X_est_fsk[i] = np.dot(F, X_est_fsk[i - 1]) + KKK * (Z[i] - A[0])
    Pk_fsk = np.dot(I - KK * H, Pk_k_1_fsk)
    if i < NN_N:
        PP_1[i] = Pk_fsk[0, 0]
        PP_2[i] = Pk_fsk[1, 1]
        PP_3[i] = Pk_fsk[2, 2]

Ve_est,Fn_est, w_est = zip(*X_est_fsk)
Ve,Fn, w = zip(*X_est)

#################Графики

figure1 = plt.figure()
plt.plot(dV)

plt.ylabel("dV, м/с")
plt.xlabel("c")

figure2 = plt.figure()
plt.plot(Fn)
plt.ylabel("Fn, рад")
plt.xlabel("c")

figure3 = plt.figure()
plt.plot(wdr)
plt.ylabel("wdr, рад/с")
plt.xlabel("c")


fig, ax1 = plt.subplots()
ax1.plot(Z1, label='График 1', color='g')
ax2 = ax1.twiny()
ax2.plot(dV, label='График 2', color='r')
# plt.plot(Zkdv1)
fig.supxlabel('Zk1, рад/с')
fig.supylabel('с')


fig1, ax1 = plt.subplots()
ax1.plot(Z2, color='g')
ax2 = ax1.twiny()
ax2.plot(dV, color='r')
fig1.supxlabel('Zk2, рад/с')
fig1.supylabel('с')

fig2, ax1 = plt.subplots()
ax1.plot(Z3, color='g')
ax2 = ax1.twiny()
ax2.plot(dV, color='r')
fig2.supxlabel('Zk3, рад/с')
fig2.supylabel('с')

fig11, ax11 = plt.subplots()
ax11.plot(Ve_est, label='График 1', color='g')
ax2 = ax11.twiny()
ax2.plot(dV, label='График 2', color='r')
fig.supxlabel('Ve')
fig.supylabel('с')

fig12, ax12 = plt.subplots()
ax12.plot(Fn_est, label='График 1', color='g')

ax2 = ax12.twiny()
ax2.plot(Fn, label='График 2', color='r')

fig.supxlabel('Ve')
fig.supylabel('с')

fig13, ax13 = plt.subplots()
ax13.plot(w_est, label='График 1', color='g')
plt.xlim(150,5000)
plt.ylim (-2e-7,2e-7)
ax2 = ax13.twiny()
ax2.plot(wdr, label='График 2', color='r')
plt.xlim((150,5000))
plt.ylim (-2e-7,2e-7)
fig.supxlabel('Ve')
fig.supylabel('с')


plt.show()
