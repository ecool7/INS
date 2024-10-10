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



# figure1 = plt.figure()
# plt.plot(dV)
#
# plt.ylabel("dV, м/с")
# plt.xlabel("c")
#
# figure2 = plt.figure()
# plt.plot(Fn)
# plt.ylabel("Fn, рад")
# plt.xlabel("c")
#
# figure3 = plt.figure()
# plt.plot(wdr)
# plt.ylabel("wdr, рад/с")
# plt.xlabel("c")
#
#
# fig, ax1 = plt.subplots()
# ax1.plot(Zkdv1, label='График 1', color='g')
# ax2 = ax1.twiny()
# ax2.plot(dV, label='График 2', color='r')
# # plt.plot(Zkdv1)
# fig.supxlabel('Zk1, рад/с')
# fig.supylabel('с')
#
#
# fig1, ax1 = plt.subplots()
# ax1.plot(Zkdv2, color='g')
# ax2 = ax1.twiny()
# ax2.plot(dV, color='r')
# fig1.supxlabel('Zk2, рад/с')
# fig1.supylabel('с')
#
#
# fig2, ax1 = plt.subplots()
# ax1.plot(Zkdv3, color='g')
# ax2 = ax1.twiny()
# ax2.plot(dV, color='r')
# fig2.supxlabel('Zk3, рад/с')
# fig2.supylabel('с')
#
# plt.show()

#
#
# #####kalman####
# T = 1e-1
# beta = 1E-11
#
# t = 10000
# N =int(t / T)
# Rz =6371000
# g = 9.80665
#
# ##Матрицы
# R = 1
# Q = 1e-15
# X0 = np.array([0,0,0])
# H = np.array([1, 0, 0])
# G = np.array([[0, 0, 0],[0, 0, 0],[0, 0, 1]])
# Xk =X0
# I = np.eye(3)
# A = np.array([[0, g, 0], [-1/Rz, 0, 1], [0, 0, 1-beta]])
# F = I + A*T + np.dot(A,A)*T**2/2
# est_X0 = np.array([[0], [0], [0]])
# est_X_prev = est_X0
# P0 = np.eye(3)
# P_prev = P0
# Nuk = np.array([[0], [0], [0]])
#
#
#
# prcnt10 = N // 100
# p_count = 0
# dV = np.zeros((N, 1))
# F = np.zeros((N, 1))
# wdr = np.zeros((N, 1))
#
# for i in range(N):
#     # est_X_trans = np.dot(F.T, est_X_prev.T)
#     Nuk = Zkdv1 - np.dot(H, np.dot(F, est_X_prev))
#     P_trans = np.dot(F, np.dot(P_prev, F.T)) + np.dot(G, np.dot(Q, G.T))
#     Kk = np.dot(P_trans, H.T) / (np.dot(H, np.dot(P_trans, H.T)) + R)
#     est_X_next =  np.dot(Kk, Nuk)
#     P_next = np.dot(np.eye(3) - np.dot(Kk, H), P_trans)
#
#     dV[i, 0] = est_X_prev[0, 0]
#     F[i, 0] = est_X_prev[1, 0]
#     wdr[i, 0] = est_X_prev[2, 0]
#
#     est_X_prev = est_X_next
#     P_prev = P_next
#
#     if i % prcnt10 == 0:
#         print(f'Filtering complete {p_count + 1}%')
#         p_count += 1
#
# time = np.arange(0, t, T)
#
# # Plot results
# plt.figure(figsize=(10, 6))
# plt.plot(time, dV, 'b', label='Estimate')
# plt.grid(True)
# plt.title('Velocity Error Estimate')
# plt.xlabel('Time (s)')
# plt.ylabel('Velocity Error (m/s)')
# plt.legend()
#
# plt.figure(figsize=(10, 6))
#
# plt.plot(time, F, 'r', label='Truth')
# plt.grid(True)
# plt.title('Orientation Error Estimate')
# plt.xlabel('Time (s)')
# plt.ylabel('Orientation Error (rad)')
# plt.legend()
#
# plt.figure(figsize=(10, 6))
#
# plt.plot(time, wdr, 'r', label='Truth')
# plt.grid(True)
# plt.title('Drift Error Estimate')
# plt.xlabel('Time (s)')
# plt.ylabel('Drift Error (rad/s)')
#
#
# plt.show()



import numpy as np
import random
import matplotlib.pyplot as plt

g = 9.81
T = 1
R = 6.3781366*10**6 # м
N = 5000
A = 0.001
W_DR_N = 4.8*10**(-8) # рад/с
beta = 0.5*9*10**(-8)
Vk = 0.1

F = np.array([[1, -g*T, 0], [T/R, 1, T], [0, 0, 1-beta*T]])
#G = [0, 0, A*T*(2*beta)**0.5]
G = np.array([0, 0, 1])#np.array([0, 0, A*T*(2*beta)**0.5])
#H = [1, 0, 0]
H = np.array([1, 0, 0])
# I = np.ones((3,3))
I = np.eye(3)
X = np.zeros((N,3)) # состояние
#X[0] = [0, 0, 0]
Z = np.zeros(N) # измерение
Z[0] = np.dot(H, X[0]) + random.random()
Q = beta**-0.5
R = Vk**2
X_est = np.zeros((N,3)) # с домиком - оценка
Pk = np.zeros((3,3)) # 3x3
Pk_k_1 = np.zeros((3,3)) # 3x3 Pk/k-1

T=np.arange(1,N)
NN_N=200
PP_1 = np.zeros(NN_N)
PP_2 = np.zeros(NN_N)
PP_3 = np.zeros(NN_N)

# Близость Домика к "Без Домика"

def delta_Bd(A:np.array, Est:np.array,N):
  Sqr_sum = 0
  for i in range(N,len(Est)):
    if Est[i] != 0:
      Sqr_sum = Sqr_sum + (A[i] - Est[i])**2
      N = N + 1
  Sqr_sum = Sqr_sum/N
  return Sqr_sum


  # Прямой калман
def StraightK():

  Q_fsk = 1e-13
  X_est_fsk = np.zeros((N,3)) # с домиком - оценка
  Pk_fsk = np.diag([10,1e-6,1e-9])# 3x3
  Pk_k_1_fsk = np.zeros((3,3)) # 3x3 Pk/k-1
  KK = np.dot(Pk_k_1_fsk, np.transpose(H)) / np.transpose(np.dot(H, np.dot(Pk_k_1_fsk, np.transpose(H))) + R) # усиление

  for i in range(1, N):
      Pk_k_1_fsk = np.dot(F, np.dot(Pk_fsk, np.transpose(F))) + np.dot(G , Q*G[:, np.newaxis])
      KK = np.dot(Pk_k_1_fsk, H[:, np.newaxis]) / (np.dot(H, np.dot(Pk_k_1_fsk, H[:, np.newaxis])) + R)
      X_old_fsk = X_est_fsk[i-1]
      A = np.dot(H, np.dot(F, X_old_fsk[:, np.newaxis]))
      KKK=np.array([KK[0,0],KK[1,0],KK[2,0]])
      X_est_fsk[i] = np.dot(F,X_est_fsk[i-1]) +  KKK*(Z[i] - A[0])
      Pk_fsk = np.dot(I-KK*H, Pk_k_1_fsk)
      if i < NN_N:
        PP_1[i] = Pk_fsk[0,0]
        PP_2[i] = Pk_fsk[1,1]
        PP_3[i] = Pk_fsk[2,2]

  return X_est_fsk


# Обратный калман
def BackwardK(X_Input:np.array):

  Q_bwk = 1e-13
  X_est_bwk = np.zeros((N,3)) # с домиком - оценка
  X_est_bwk[0] = X_Input[N-1]
  Pk_bwk = np.diag([10,1e-6,1e-9])# 3x3
  Pk_k_1_bwk = np.zeros((3,3)) # 3x3 Pk/k-1
  KK = np.dot(Pk_k_1_bwk, np.transpose(H)) / np.transpose(np.dot(H, np.dot(Pk_k_1_bwk, np.transpose(H))) + R) # усиление

  F_bwk = np.linalg.inv(F)

  for i in range(1, N):
      Pk_k_1_bwk = np.dot(F_bwk, np.dot(Pk_bwk, np.transpose(F_bwk))) + np.dot(G , Q*G[:, np.newaxis])
      KK = np.dot(Pk_k_1_bwk, H[:, np.newaxis]) / (np.dot(H, np.dot(Pk_k_1_bwk, H[:, np.newaxis])) + R)
      X_old_bwk = X_est_bwk[i-1]
      A = np.dot(H, np.dot(F_bwk, X_old_bwk[:, np.newaxis]))
      KKK=np.array([KK[0,0],KK[1,0],KK[2,0]])
      X_est_bwk[i] = np.dot(F_bwk,X_est_bwk[i-1]) +  KKK*(Z[N-i] - A[0])
      Pk_bwk = np.dot(I-KK*H, Pk_k_1_bwk)
      if i < NN_N:
        PP_1[i] = Pk_bwk[0,0]
        PP_2[i] = Pk_bwk[1,1]
        PP_3[i] = Pk_bwk[2,2]

  return X_est_bwk

#tail cut
def Tail_f(PP:np.array):
  Win = 20
  eps = 0.1
  Mean = 0
  SWin = 0
  for i in range(len(PP)):
    if i < Win:
      SWin = SWin+PP[i]
    if i >= Win:
      Mean = SWin/20
      SWin = SWin-PP[i-Win]+PP[i]
      if abs((PP[i]-Mean)/Mean)<=eps:
        return i
        break

for k in range(1, N):
  X[k] = np.dot(F, X[k-1]) + np.dot(G, (-1)**k*random.random())
  Z[k] = np.dot(H, X[k]) + (-1)**k*Vk*random.random()


# Построение графиков

Ve,Fn, w = zip(*X)

x_label = "Время, сек"

fig1, (fig1_ax1, fig1_ax2, fig1_ax3, fig1_ax4) = plt.subplots(nrows=4, ncols=1, figsize=(18, 10))

fig1_ax1.plot(Ve,color = 'k')
fig1_ax1.set_ylabel("Ошибка по скорости, м/с")
fig1_ax1.set_xlabel(x_label)
fig1_ax1.grid(which='major', linewidth=2)
fig1_ax1.grid(which='minor')

fig1_ax2.plot(Fn, color = 'k')
fig1_ax2.set_ylabel("Фn, рад")
fig1_ax2.set_xlabel(x_label)
fig1_ax2.grid(which='major', linewidth=2)
fig1_ax2.grid(which='minor')

fig1_ax3.plot(w, color = 'k')
fig1_ax3.set_ylabel("Дрейф, рад/с")
fig1_ax3.set_xlabel(x_label)
fig1_ax3.grid(which='major', linewidth=2)
fig1_ax3.grid(which='minor')

fig1_ax4.plot(Z, color = 'k')
fig1_ax4.set_ylabel("Измерение, м/с")
fig1_ax4.set_xlabel(x_label)
fig1_ax4.grid(which='major', linewidth=2)
fig1_ax4.grid(which='minor')

plt.show()