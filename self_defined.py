import numpy as np
from gmpy2 import is_prime
import gmpy2


def xxrange(num, row):
    l = []
    for i in range(num):
        if i != row:
            l.append(i)
    return l


def Adjugate_matrix(sub_B):
    '''
    求伴随矩阵
    :param sub_B: a 2D square matrix: [[], [], ...]
    :return:
    '''
    row_num, col_num = len(sub_B), len(sub_B[0])

    # 初始化
    A_sub_B = []
    for i in range(row_num):
        A_sub_B.append([])

    for row in range(row_num):
        for col in range(col_num):
            p = pow(-1, row+1+col+1)

            kept_rows = xxrange(row_num, row)
            kept_cols = xxrange(col_num, col)
            tmp_sub_B = []
            for kept_row in kept_rows:
                c = []
                for kept_col in kept_cols:
                    c.append(sub_B[kept_row][kept_col])
                tmp_sub_B.append(c)

            tmp_D = np.linalg.det(np.array(tmp_sub_B))

            A_sub_B[row].append(gmpy2.mul(gmpy2.mpz(round(tmp_D)), p))

    return A_sub_B


def matrix_element_mod(sub_B, m):
    '''
    矩阵每个元素模m
    :param sub_B: a 2D square matrix: [[], [], ...]
    :param m: 模数
    :return: a 2D square matrix: [[], [], ...]
    '''
    row_num, col_num = len(sub_B), len(sub_B[0])

    # 初始化
    sub_B_mod_m = []
    for row in range(row_num):
        sub_B_mod_m.append([])

    for row in range(row_num):
        for col in range(col_num):
            sub_B_mod_m[row].append(gmpy2.f_mod(gmpy2.mpz(sub_B[row][col]), gmpy2.mpz(m)))
    return sub_B_mod_m


def matrix_inverse(sub_B, n):
    '''
    求矩阵模n的逆（已知sub_B一定有逆）
    :param sub_B: 一个2D矩阵: [[], [], ...]
    :param n: 一个正整数: gmpy2.mpz()
    :return: 矩阵模n的逆
    '''
    # print(np.linalg.det(sub_B))
    # det会有精度问题（矩阵元素都是int但是行列式是float, 四舍五入问题解决）

    D = gmpy2.mpz(round(np.linalg.det(np.array(sub_B))))
    D = pow(D, 1, n)
    # print('Det of sub_B: ', D)

    D_inverse = gmpy2.invert(D, n)
    # print('sub_B的行列式模{}的逆: {}'.format(n, D_inverse))

    A_star = Adjugate_matrix(sub_B)
    # print('sub_B的伴随矩阵：', A_star)

    A_star_times_D_inverse = []
    for row in A_star:
        c = []
        for row_col in row:
            c.append(gmpy2.mul(row_col, D_inverse))
        A_star_times_D_inverse.append(c)

    A_star_times_D_inverse_mod_n = matrix_element_mod(A_star_times_D_inverse, n)

    # print('行列式的逆 乘以 sub_B的伴随矩阵 模{}：{}'.format(n, A_star_times_D_inverse_mod_n))

    return np.transpose(np.array(A_star_times_D_inverse_mod_n))


def dot(A, B):
    '''
    元素都是大整数的向量点乘，防止出现类似“OverflowError: cannot fit 'mpz' into an index-sized integer”的错误
    :param A: list
    :param B: list
    :return: list
    '''
    if len(A) != len(B):
        return -1
    else:
        l = len(A)
    dot_product = 0
    for i in range(l):
        dot_product = gmpy2.add(dot_product, gmpy2.mul(A[i], B[i]))

    return dot_product

if __name__ == '__main__':
    # print(dot([1, 1, 4], [5, 6, 7]))
    x = gmpy2.t_mod_2exp(39, 2)
    print(x)
def isPrime(num):
    '''
    Determine if a number is a prime
    :param num: int
    :return: BOOL
    '''
    if num <= 1:
        return False

    for i in range(2, int(num ** 0.5) + 1):
        if num % i == 0:
            return False

    return True


def gcd(a, b):
    if a < b:
        a, b = b, a
    # print(a,b)
    if a % b == 0:
        return b
    else:
        return gcd(b, a % b)


def lcm(a, b):
    '''
    :param a: int
    :param b: int
    :return: int, least common multiple of p and q
    '''
    k = gcd(a, b)
    if k == 1:
        # print('The inputs are primes! ')
        return a*b
    else:
        return int((a*b)/k)


def egcd(a, b):
    '''
    extended gcd
    :param a: int
    :param b: int
    :return: x, y, where x*a + y*b = gcd(a, b)
    '''
    if b == 0:
        return 1, 0
    else:
        k = a // b
        r = a % b
        x1, y1 = egcd(b, r)
        x, y  = y1, x1 - k*y1
    return x, y


def multiplicative_inverse(a, b):
    '''
    求解a mod b 的逆元x，即ax = 1 mod b
    :param a: int
    :param b: int
    :return: x
    '''
    if gcd(a, b) == 1:
        # 只有互素才有逆元
        x, y = egcd(a, b)
        return x % b
    else:
        print('{} 和 {}不互素，没有逆元'.format(a, b))
        return -1


def quick_mod_pow(a, b, c):
    A = 1
    T = a % c
    while(b != 0):
        # 逻辑与；例如，1010 & 0001 = 0
        # 作用：提取出最右边的二进制位
        if (b & 1):
            A = (A*T) % c
        # 右移一位；例如，1001>>=1 => 110
        b >>= 1
        T = (T*T) % c
    # print(A)
    return A


def isInvertible(A):
    '''
    Whether A is invertible.
    :param A: a matrix
    :return: BOOL
    '''
    if A.shape[0] != A.shape[1]:
        return False
    else:
        if round(np.linalg.det(A)) == 0:
            return False
        else:
            return True
