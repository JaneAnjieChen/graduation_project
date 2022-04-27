import time
import numpy as np
from Asymmetric_Cipher import calc_B_comb_inverse, SKeyGen, SKeySplit
# 3. 子阵模逆的运算
TIMES = 10
pk, sk = SKeyGen(bit_length=1024)
partial_strong_private_keys, B = SKeySplit(SK=sk)

csp_indices = [1, 2, 3, 4, 5, 7, 8]
average_time = []

for i in range(TIMES):
    start = time.clock()
    B_comb_inverse = calc_B_comb_inverse(SK=sk, csp_indices=csp_indices, B=B)
    end = time.clock()
    average_time.append(end-start)
print(np.average(average_time))
