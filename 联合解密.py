import time
import numpy as np
from Asymmetric_Cipher import SKeyGen, SKeySplit, Enc, WKeygen, PSDec, calc_B_comb_inverse, Comb

# 7. 联合解密
TIMES = 10
pk, sk = SKeyGen(bit_length=1024)
partial_strong_private_keys, B = SKeySplit(SK=sk)
wpk, wsk = WKeygen(PK=pk)
C, W = Enc(wpk, m=51)

csp_indices = [1, 2, 3, 4, 5, 7, 8]
csp_indices = sorted(csp_indices)
shares = {}
for i in csp_indices:
    c_i = PSDec(C, partial_strong_private_keys[i])
    shares[i] = c_i
# print('shares: ', shares)

B_comb_inverse = calc_B_comb_inverse(SK=sk, csp_indices=csp_indices, B=B)

avarege_time = []
for t in range(TIMES):
    start = time.clock()
    Comb(shares, B_comb_inverse)
    end = time.clock()
    avarege_time.append(end-start)
print(np.average(avarege_time))
