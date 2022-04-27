# from self_defined import isPrime, multiplicative_inverse, quick_mod_pow
import random
import math
from Asymmetric_Cipher import prime_generator
import gmpy2

# # this is pre-set
# def get_delta(psai):
#     # (2 ** mi)<= psai <= (2 ** (mi+1))
#     mi = -1
#     while psai != 0:
#         psai >>= 1
#         mi += 1
#     print(mi)
#
#     # delta1 + delta2 <= mi
#     sum = random.randint(10, mi)
#     print(sum)
#     delta1 = random.randint(5, int(sum/2))
#     delta2 = sum - delta1
#     return delta1, delta2

def key_generation():
    import global_variable as gl
    key = random.randint(1, gmpy2.sub(gl.psai, 1))
    return key

def E(m, key):
    import global_variable as gl
    psai = gl.psai
    delta1 = gl.delta1

    pow_2_delta1 = pow(2, delta1)
    e = random.randint(0, gmpy2.sub(pow(2, gl.k), 1))
    # print('e: ', e)

    part_1 = gmpy2.t_mod(gmpy2.mul(e, pow_2_delta1), psai)
    # print('part1: ', part_1)
    part_2 = gmpy2.t_mod(gmpy2.add(m, part_1), psai)
    # print('part2: ', part_2)
    c = gmpy2.t_mod(gmpy2.mul(key, part_2), psai)
    # print(c)

    return c


def D(c, k):
    import global_variable as gl
    psai = gl.psai
    delta1 = gl.delta1

    # pow_2_delta1 = pow(2, delta1)
    inverse_k = gmpy2.invert(k, psai)
    # print('inverse_k: ', inverse_k)
    part_1 = gmpy2.t_mod(gmpy2.mul(inverse_k, c), psai)
    m = gmpy2.t_mod_2exp(part_1, delta1)
    # print(m)

    return m


if __name__ == '__main__':
    psai = prime_generator(requested_number=1, requested_length=600)[0]
    print("psai: ", psai)
    import global_variable as gl
    gl.psai = psai

    delta1, delta2 = 70, 530
    # print("delta1: ", delta1)
    # print("delta2: ", delta2)
    import time

    start = time.clock()
    k1 = key_generation()
    end = time.clock()
    print(end-start)
#     print('key: ', k1)
#
#
#     m = 19400828
#     print('message: ', m)
#
#     c = E(m, k1)
#
#     D(c, k1)
