import time
import numpy as np
from Asymmetric_Cipher import SKeyGen, SKeySplit
# 2. 强私钥的分发
TIMES = 10

pk, sk = SKeyGen(bit_length=1024)

import global_variable as gl
print(gl.t)

average_time = []
for t in range(TIMES):
    start = time.clock()
    partial_strong_private_keys, B = SKeySplit(SK=sk)
    end = time.clock()
    average_time.append(end-start)

print(np.average(average_time))
