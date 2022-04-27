# 先定义函数模块，再用类（角色）调用这些模块
import random
import math
import numpy as np
import itertools
from self_defined import matrix_inverse
# from gmpy2 import is_prime, mpz, next_prime, add, t_div
import gmpy2


def prime_generator(requested_number, requested_length):
    '''
    :param requested_number: the number of primes requested
    :param requested_length: the requested bit length of the primes
    :return: randomly selected primes, where (prime-1)/2 is also prime
    '''
    primes = []
    start = pow(gmpy2.mpz(2), gmpy2.mpz(requested_length))

    tmp = gmpy2.next_prime(start)
    while len(primes) < requested_number:
        tmp1 = gmpy2.t_div(gmpy2.sub(tmp, gmpy2.mpz(1)), gmpy2.mpz(2))
        if gmpy2.is_prime(tmp1) and tmp1 not in primes:
            primes.append(tmp)
            tmp = gmpy2.next_prime(tmp)
        else:
            tmp = gmpy2.next_prime(tmp)

    return primes


def SKeyGen(bit_length=1024):
    '''
    :param bit_length: the bit length
    :return: the common pubilc key pk = (N, g); the strong private key sk = lamda
    '''

    # Select 2 large enough primes, where (prime-1)/2 is also prime
    # the primes is of "bit_length" bits
    [p, q] = prime_generator(requested_number=2, requested_length=bit_length)
    print("p {}".format(p))
    print("q {}".format(q))

    # N
    N = gmpy2.mul(p, q)
    # print('N: ', N)

    p1 = gmpy2.t_div(gmpy2.sub(p, gmpy2.mpz(1)), gmpy2.mpz(2))
    q1 = gmpy2.t_div(gmpy2.sub(q, gmpy2.mpz(1)), gmpy2.mpz(2))

    # lamda
    lamda = gmpy2.mul(p1, q1)
    # print('lamda: {}'.format(lamda))

    N2 = gmpy2.mul(N, N)

    # g不可以取1，因为1和任何数都不互素

    # 证明可得g可以这么取
    # g = pow(2, gmpy2.mul(2, N), N2)
    a = 3
    import global_variable as gl
    gl.N = N
    gl.a = a
    gl.N2 = N2

    return (N, a), lamda


def SKeySplit(SK):
    '''
    Split the strong privte key into n parts.
    :return: a list of private keys derived from the strong private key
    '''
    # a1 = lamda * (lamda ** (-1) mod (N*N))
    import global_variable as gl
    n = gl.n
    t = gl.t
    N = gl.N
    lamda = SK
    N2 = gl.N2

    MAX_INT = 100
    if MAX_INT < N:
        B_HIGH = MAX_INT

    flag = 1
    while flag == 1:
        flag = 0
        # 使得B的所有元素为正整数
        B = np.random.randint(1, B_HIGH, (n, t))
        # n行中选t行
        for rows in itertools.combinations(range(0, n), t):
            rows = list(rows)
            sub_B = B[rows, :]
            # sub_B的行列式
            D = gmpy2.mpz(round(np.linalg.det(sub_B)))
            D_mod = gmpy2.t_mod(D, N2)
            # 是否在模上可逆
            if gmpy2.gcd(D_mod, N) != 1 or gmpy2.gcd(D_mod, lamda) != 1:
                flag = 1
                break
            # if gmpy2.gcd(D_mod, 2) != 1:
            #     flag = 1
            #     break

    # B需要保密
    # import global_variable as gl
    # gl.B = B
    # print('Matrix B: {}'.format(B))

    if n < t:
        return -1

    lamda_inverse = gmpy2.invert(lamda, N2)
    # print("lamda's inverse: ", lamda_inverse)

    a1 = gmpy2.mul(lamda, lamda_inverse)
    vector_a = [a1]

    # randomly selected {a2, a3, ..., at}
    for i in range(t-1):
        vector_a.append(random.randint(0, B_HIGH))
    # vector_a = np.array(vector_a)
    # print('vector_a: ', vector_a)

    vector_lamda = np.dot(B, vector_a)
    # print(type(vector_lamda))

    vector_lamda_list = vector_lamda.tolist()
    # print(type(vector_lamda_list))

    # return the partial strong private keys
    return vector_lamda_list, B


def calc_B_comb_inverse(SK, csp_indices, B):
    '''
    可信计算方计算B_comb的逆
    :param csp_indices: 参与comb计算的CSP的索引列表
    :param B: 随机生成的B
    :return: B_comb_inverse
    '''
    import global_variable as gl

    N = gl.N
    lamda = SK
    # print('B: ', B)
    B_comb = B[csp_indices, :]
    # print('B_comb: ', B_comb)
    B_comb = B_comb.tolist()

    B_comb_inverse = matrix_inverse(B_comb, gmpy2.mul(lamda, N))

    return B_comb_inverse


def WKeygen(PK):
    '''
    Weak key generation
    :param PK: strong public key
    :return: weak public key, weak private key
    '''
    import global_variable as gl
    N2 = gl.N2
    (N, a) = PK
    omiga = gmpy2.mpz(random.randint(0, gmpy2.t_div(N, 4)))
    g = pow(a, gmpy2.mul(2, N), N2)
    h = pow(g, omiga, N2)
    return (N, a, h), omiga


def Enc(pkj, m):
    '''
    Encryption with a week public key
    :param pkj: weak public key pkj
    :param m: message, m < N
    :return: (C, W)
    '''
    (N, a, h) = pkj
    import global_variable as gl
    N2 = gl.N2

    r = gmpy2.mpz(random.randint(0, gmpy2.t_div(N, 4)))
    # print('r: ', r)

    C1 = pow(h, r, N2)
    # print(C1)
    C2 = pow(gmpy2.add(1, gmpy2.mul(m, N)), 1, N2)
    # print(C2)
    # print(C1*C2)
    C = pow(gmpy2.mul(C1,C2), 1, N2)
    g = pow(a, gmpy2.mul(2, N), N2)
    W = pow(g, r, N2)

    return (C, W)


def WDec(C, W, wskj):
    '''
    Decryption with a weak private key
    :param C: 含有真正信息的密文
    :param W: 辅助弱私钥解密的部分
    :param wskj: 弱私钥
    :return: message m
    '''
    omiga = wskj
    import global_variable as gl
    N = gl.N
    N2 = gl.N2

    # 计算逆而不是除法
    C_part1 = pow(C, 1, N2)
    # print('分母：', C_part1)
    C_part2 = gmpy2.invert(pow(W, omiga, N2), N2)
    # print('分子：', C_part2)

    x = pow(gmpy2.mul(C_part1, C_part2), 1, N2)
    # print(x)

    return gmpy2.t_div(gmpy2.sub(x, 1), N)


def SDec(C, SK):
    '''
    Decryption with a strong private key.
    :param C: 含有真正信息的密文
    :param SK: strong private key
    :return: message m
    '''
    import global_variable as gl
    N = gl.N
    lamda = SK
    N2 = gl.N2
    x = pow(C, lamda, N2)
    # print(x)

    L = gmpy2.t_div(gmpy2.sub(x, 1), N)
    # print('L: ', L)

    lamda_inverse = gmpy2.invert(lamda, N2)
    # print('lamda_inverse: ', lamda_inverse)

    return pow(gmpy2.mul(L, lamda_inverse), 1, N)


def PSDec(C, partial_SK):
    '''
    Partial decryption with a partial strong private ket
    :param C: 含有真正信息的密文
    :param partial_SK: 经过SKeySplit的强私钥
    :return: a decryption share c_i
    '''
    import global_variable as gl
    N = gl.N
    lamda_i = partial_SK
    N2 = gl.N2

    if lamda_i >= 0:
        c_i = pow(C, lamda_i, N2)
    else:
        c_i = pow(gmpy2.invert(C, N2), (-lamda_i), N2)
    # print('decryption share:', c_i)
    return c_i



def Comb(shares, B_comb_inverse):
    '''
    :param shares: {}
    :param B_comb_inverse: B_comb的模逆
    :return:
    '''

    import global_variable as gl
    N = gl.N
    N2 = gl.N2

    if len(shares.keys()) < gl.t:
        print('NOT ENOUGH SHARES.')
        return -1

    # print('B_comb_inverse: ', B_comb_inverse)

    # 先行后列 B_comb_inverse[0][]
    product = 1
    B_comb_inverse_row0_list = B_comb_inverse[0].tolist()
    # print(B_comb_inverse_row0_list)
    # print('~'*100)

    count = 0
    for i in shares.keys():
        b = B_comb_inverse_row0_list[count]
        # 大于0，不计算逆，直接指数运算
        if b >= 0:
            tmp = pow(shares[i], b, N2)
        # 小于0，计算逆，再指数运算
        else:
            tmp = pow(gmpy2.invert(shares[i], N2), (-b), N2)
        # print(tmp)
        product = gmpy2.mul(product, tmp)
        count += 1
    # print('product: ', product)
    x = pow(product, 1, N2)
    # x = product % N2
    # print('x: ', x)

    L = gmpy2.t_div(gmpy2.sub(x, 1), N)

    return int(L)


if __name__ == '__main__':
    # import time
    # start = time.clock()
    # print(start)

    strong_public_key, strong_private_key = SKeyGen(bit_length=1024)
    print('PK: {}, SK: {}'.format(strong_public_key, strong_private_key))

    # end0 = time.clock()
    # print(end0 -start)
    # start = end0

    # partial_strong_private_keys, B = SKeySplit(SK=strong_private_key)
    # print('partial_strong_private_keys: ', partial_strong_private_keys)
    #
    # # end0 = time.clock()
    # # print(end0 -start)
    # # start = end0
    #
    # weak_public_key, weak_private_key = WKeygen(PK=strong_public_key)
    # print('pkj: {}, wskj: {}'.format(weak_public_key, weak_private_key))
    #
    # # end0 = time.clock()
    # # print(end0 - start)
    # # start = end0
    #
    # # # 攻击者
    # # N = weak_public_key[0]
    # # a = weak_public_key[1]
    # # h = weak_public_key[2]
    # # attack_get_wpk = gmpy2.log(h)/gmpy2.log(gmpy2(a, gmpy2.mul(2, N)))
    # # print(attack_get_wpk)
    #
    # m = prime_generator(requested_number=1, requested_length=600)[0]
    # C, W = Enc(weak_public_key, m)
    # print('C: {}, W: {}'.format(C, W))
    # print('-'*100)
    #
    # m1 = WDec(C, W, weak_private_key)
    # print(m1)
    #
    # m2 = SDec(C, SK=strong_private_key)
    # print(m2)
    #
    # import global_variable as gl
    # shares = {}
    # csp_indices = [1, 2, 3, 4, 5, 7, 8]
    # csp_indices = sorted(csp_indices)
    # for i in csp_indices:
    #     c_i = PSDec(C, partial_strong_private_keys[i])
    #     shares[i] = c_i
    # print('shares: ', shares)
    #
    # # shares = dict([(csp_index, shares[csp_index] for csp_index in csp_indices)])
    #
    # B_comb_inverse = calc_B_comb_inverse(SK=strong_private_key, csp_indices=csp_indices, B=B)
    #
    # m3 = Comb(shares, B_comb_inverse)
    # print(m3)
