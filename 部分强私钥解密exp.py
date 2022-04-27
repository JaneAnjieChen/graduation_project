import time
import numpy as np
from Asymmetric_Cipher import SKeyGen, SKeySplit, Enc, WKeygen, PSDec
# # 6. 部分强私钥的解密
# TIMES = 10
# pk, sk = SKeyGen(bit_length=1024)
# partial_strong_private_keys, B = SKeySplit(SK=sk)
# wpk, wsk = WKeygen(PK=pk)
# C, W = Enc(wpk, m=51)
#
# diff_i_time = []
# for i in range(len(partial_strong_private_keys)):
#     average_time = []
#
#     for t in range(TIMES):
#         start = time.clock()
#         c_i = PSDec(C, partial_strong_private_keys[i])
#         end = time.clock()
#         average_time.append(end-start)
#
#     diff_i_time.append(np.average(average_time))
#
# print(diff_i_time)

diff_i_time = [0.03320628999999999, 0.03280344, 0.03301261999999999, 0.032637719999999995, 0.032644530000000005, 0.03267603000000001, 0.03333572999999999, 0.032887629999999966, 0.032567739999999956, 0.03288868000000007]
indices = [0,1,2,3,4,5,6,7,8,9]
import matplotlib.pyplot as plt
from scipy import optimize

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

plt.plot(indices, diff_i_time, '+', color='b')

# 拟合
def f_1(x, A, B):
    return A * x + B
A1, B1 = optimize.curve_fit(f_1, indices, diff_i_time)[0]
print(A1, B1)
x1 = np.arange(indices[0], indices[len(indices)-1], 0.01)
y1 = round(A1, 4)*x1 + B1
plt.plot(x1, y1, color='r', label='拟合')
plt.yticks(np.linspace(0, 0.1, 2))
plt.xticks(np.linspace(0, 9, 10))
plt.xlabel('序号')
plt.ylabel('用时（秒）')

plt.legend()
plt.show()
