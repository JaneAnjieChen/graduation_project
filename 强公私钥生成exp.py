import time
import numpy as np
# 1. 强公钥的生成（模块）
from Asymmetric_Cipher import SKeyGen

bit_lengths = [100, 200, 300, 400, 500, 600, 700, 800, 1000, 1100]
average_time = []
TIMES = 10
for bit_length in bit_lengths:
    elapses = []
    for t in range(TIMES):
        start = time.clock()
        pk, sk = SKeyGen(bit_length=bit_length)
        end = time.clock()
        elapses.append(end - start)
    average_time.append(np.mean(elapses))

print(average_time)

# import matplotlib.pyplot as plt
# plt.rcParams['font.sans-serif']=['SimHei']
# plt.rcParams['axes.unicode_minus'] = False
# plt.bar(bit_lengths, average_time, width=1)
# x_ticks = ['256', '512', '1024', '2048']
# plt.xticks()
# plt.xlabel('素数比特长度')
# plt.ylabel('平均运行时间/秒')
# plt.show()
