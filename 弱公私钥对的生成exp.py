import time
import numpy as np
from Asymmetric_Cipher import WKeygen, SKeyGen, SKeySplit
# 4. 弱公私钥对的生成
TIMES = 10
pk, sk = SKeyGen(bit_length=1024)

average_time = []

for i in range(TIMES):
    start = time.clock()
    wpk, wsk = WKeygen(PK=pk)
    end = time.clock()
    average_time.append(end-start)
print(np.average(average_time))
