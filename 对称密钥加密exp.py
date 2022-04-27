import time
import numpy as np
from Asymmetric_Cipher import prime_generator
from Symmetric_Cipher import key_generation, E

# 1. 相同参数，相同k，加密不同长度的信息
psai = prime_generator(requested_number=1, requested_length=600)[0]

import global_variable as gl
gl.psai = psai
gl.delta1 = 70
gl.delta2 = 530

key = key_generation()

TIMES = 10

bit_length = []
a = []
for b in range(0, 35):
    average_time = []

    m = pow(2, b)

    bit_length.append(b+1)

    for i in range(TIMES):
        start = time.clock()
        E(m=m, key=key)
        end = time.clock()
        average_time.append(end-start)

    a.append(np.average(average_time))

print(a)



import matplotlib.pyplot as plt
from scipy import optimize

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

plt.plot(bit_length, a, '+', color='b')

# 拟合
def f_1(x, A, B):
    return A * x + B
A1, B1 = optimize.curve_fit(f_1, bit_length, a)[0]
print(A1, B1)
x1 = np.arange(bit_length[0], bit_length[len(bit_length)-1], 0.01)
y1 = round(A1, 4)*x1 + B1
plt.plot(x1, y1, color='r', label='拟合')

plt.yticks(np.linspace(0, 0.01, 2))

plt.legend()
plt.xlabel('信息的比特长度')
plt.ylabel('加密时间（毫秒）')
plt.show()
