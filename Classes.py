from Asymmetric_Cipher import SKeyGen, SKeySplit, WKeygen, Enc, PSDec, Comb, calc_B_comb_inverse
from self_defined import dot
from Symmetric_Cipher import prime_generator, key_generation, E
import math
import time
import random
import gmpy2
import numpy as np
import itertools
import global_variable as gl

k = gl.k

def string_to_vector(feature_string):
    # L = 16
    # 元素范围0-3
    feature_vector = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0]
    return feature_vector

class HMC:
    def __init__(self):
        pk, sk = SKeyGen(bit_length=1024)
        self.pk = pk
        self.sk = sk

    def initialize_SSHE_para(self):
        psai = prime_generator(requested_number=1, requested_length=600)[0]

        import global_variable as gl
        gl.psai = psai
        # 向量元素范围
        gl.ita = 2
        # 向量维度
        gl.L = 16

        # 根据条件式推到得到
        # 所以结合以上两个不等式：
        delta1 = 70
        delta2 = 530
        # Make sure: delta1 + delta2 < log(psai, 2))
        if pow(2, delta1 + delta2) < psai:
            gl.delta1 = delta1
            gl.delta2 = delta2
        else:
            print('SSHE paremeters initialization failed.')

    def HMC_generate_SSHE_key(self):
        import global_variable as gl
        key = key_generation()
        self.key = key

        wpk, wsk = WKeygen(PK=self.pk)
        self.wpk = wpk
        self.wsk = wsk
        (C, W) = Enc(pkj=wpk, m=key)
        self.C = C
        self.W = W

    def publish(self):
        self.initialize_SSHE_para()
        self.HMC_generate_SSHE_key()
        return self.wpk, (self.C, self.W)

    def HMC_SkeySplit(self):
        psk_list, B = SKeySplit(SK=self.sk)
        self.psk_list = psk_list
        self.B = B

    def HMC_calc_B_comb_inverse(self):
        self.B_comb_inverts = {}
        for rows in itertools.combinations(range(0, gl.n), gl.t):
            rows = list(rows)
            B_comb_invert = calc_B_comb_inverse(self.sk, rows, self.B)
            self.B_comb_inverts[tuple(rows)] = B_comb_invert

    def Doctor_verification(self, feature_string):
        # 如果验证通过
        return True

    def feature_string_to_vector(self, feature_string):
        # 从字符转化为32维、元素范围为{0, ..., 2^32 - 1}的向量
        feature_vector = string_to_vector(feature_string)
        return feature_vector

    def HMC_SSHE_encrypt_feature_vector(self, feature_vector_list):
        encrypted_vector_list = []
        vector_sq_sum = 0

        for e in feature_vector_list:
            encrypted_e = E(m=e, key=self.key)
            encrypted_vector_list.append(encrypted_e)
            vector_sq_sum = gmpy2.add(vector_sq_sum, gmpy2.mul(e, e))
        return encrypted_vector_list, vector_sq_sum


class DRP:
    def __init__(self):
        import global_variable as gl
        pk = (gl.N, gl.a)

        wpk, wsk = WKeygen(PK=pk)
        self.wpk = wpk
        self.wsk = wsk
        self.doctor_info = {}
        self.patient_info = {}
        self.k_dec1 = None

    def publish(self):
        return self.wpk

    def storage(self, info_index, info, IsDoctorInfo=False, IsPatientInfo=False, k_dec1=False):
        if IsDoctorInfo:
            self.doctor_info[info_index] = info
        elif IsPatientInfo:
            self.patient_info[info_index] = info
        elif k_dec1:
            self.k_dec1 = info

    def calc_r_patient_key(self, patient_index):
        import global_variable as gl
        self.r = random.randint(0, gl.k)

        C_k_patient, W_k_patient = self.patient_info[patient_index][2]

        N2 = gl.N2
        self.r_C_k_patient = pow(C_k_patient, self.r, N2)
        self.W_k_patient = W_k_patient

    def calc_sim_values_and_recommend_doctors(self, patient_index):
        import global_variable as gl
        # r_inverse = gmpy2.invert(self.r, gl.psai)
        k_dec = gmpy2.t_mod(gmpy2.mul(self.k_dec1, self.r), gl.psai)
        self.k_dec = k_dec
        A_ = np.array(self.patient_info[patient_index][0])
        sim_ = self.patient_info[patient_index][1]
        self.recommend_doctors = []

        for doctor_index, (doctor_vector_, doctor_sq_sum) in self.doctor_info.items():
            # start = time.clock()
            # print(start)
            I_ = dot(A_, doctor_vector_)
            # print('密文点乘：', I_)
            I = gmpy2.t_mod_2exp(gmpy2.t_mod(gmpy2.mul(k_dec, I_), gl.psai), gl.delta1)
            # print('由密文点乘计算出的明文点乘: ', I)
            fi = sim_ - doctor_sq_sum + gmpy2.mul(2, I)
            # print('fi: ', fi)
            if fi >= 0:
                self.recommend_doctors.append(doctor_index)
            # end = time.clock()
            # print(end)
            # print("多一个医生计算用时：{}".format(end-start))

    def calc_r_product_to_csp(self, product, count, wpk_hmc):
        import global_variable as gl
        r1 = gmpy2.mpz(random.randint(0, gmpy2.t_div(gl.N, 2)))
        r1_encrypted = Enc(pkj=wpk_hmc, m=r1)[0] # 只要C
        self.r1 = r1
        self.v_count = count
        return gmpy2.t_mod(gmpy2.mul(product, r1_encrypted), gl.N2)

    def calc_doctor_value(self, m):
        return gmpy2.div(gmpy2.sub(m[0], self.r1), self.v_count)


class EHR:
    def __init__(self):
        self.record = {} # format {patient_index:{doctor_index_1: [s1, s2, s3....], doctor_index_2: [...]}, ...}

    def storage(self, doctor_index, patient_index, v):
        # v format: (v1, v_hmc)
        if patient_index in self.record:
            if doctor_index in self.record[patient_index]:
                self.record[patient_index].append(v)
            else:
                self.record[patient_index] = [v]
        else:
            self.record[patient_index] = {}
            self.record[patient_index][doctor_index] = [v]

    # 部分模拟
    def extract_v_from_record(self):
        doctor_v = {} # format {doctor_index: [( (s_C,s_W),(s_C,s_W) )...]}
        for patient_index, value in self.record.items():
            for doctor_index, scores in value.items():
                if doctor_index in doctor_v:
                    doctor_v[doctor_index] += scores
                else:
                    doctor_v[doctor_index] = scores
        self.doctor_v = doctor_v

    def mul_v(self, doctor_index):
        import global_variable as gl
        product = 1
        count = 0
        for i in self.doctor_v[doctor_index]:
            # 只取hmc加密结果的C
            product = gmpy2.t_mod(gmpy2.mul(product, i[1][0]), gl.N2)
            count += 1
        return product, count


class CSP:
    def __init__(self, index, psk):
        self.index = index
        self.psk = psk
        self.to_be_partially_decrypted = [] # [(C, W),()...]

    def storage(self, info):
        self.to_be_partially_decrypted += info

    def partial_decryption(self):
        self.shares = []
        for t in self.to_be_partially_decrypted:
            try:
                if len(t) > 1:
                    C = t[0]
            except:
                C = t
            share = PSDec(C=C, partial_SK=self.psk)
            self.shares.append(share)

        self.to_be_partially_decrypted = []


class chosen_CSP:
    def __init__(self, index, B_comb_inverts):
        self.index = index
        self.shares_to_be_combined = {} # {index: (int)} / {index: (int, int)}
        self.B_comb_inverts = B_comb_inverts

    def storage(self, info, IsShares=False):
        if IsShares:
            self.shares_to_be_combined.update(info)

    def combine_shares(self, csp_indices):
        # one CSP
        for csp_index in self.shares_to_be_combined.keys():
            one_csp = csp_index
            break

        # 找到模逆
        B_comb_invert = self.B_comb_inverts[tuple(csp_indices)]

        l = len(self.shares_to_be_combined[one_csp])
        self.result = []
        for i in range(l):
            shares = {}
            for csp_index, share in self.shares_to_be_combined.items():
                shares[csp_index] = share[i]
            self.result.append(Comb(shares, B_comb_invert))

        self.shares_to_be_combined = {}

    def calc_result_to_drp(self):
        # len(result) == 2
        import global_variable as gl
        if len(self.result) == 2:
            product = gmpy2.mul(gmpy2.mpz(self.result[0]), gmpy2.mpz(self.result[1]))
            product_inverse = gmpy2.invert(product, gl.psai)
            self.k_dec1 = E(m=1, key=product_inverse)


class Doctor:
    def __init__(self, feature_string, index=None, feature_vector=None):
        self.feature_string = feature_string
        self.index = index
        self.feature_vector = feature_vector


class Patient:
    def __init__(self, index, feature_string=None):
        import global_variable as gl
        pk = (gl.N, gl.a)

        wpk, wsk = WKeygen(PK=pk)
        self.wpk = wpk
        self.wsk = wsk
        self.index = index
        self.feature_string = feature_string

    def publish(self):
        return self.wpk

    def feature_string_to_vector(self):
        self.feature_vector = string_to_vector(self.feature_string)

    def Patient_generate_SSHE_key(self):
        key = key_generation()
        self.key = key

        (C, W) = Enc(pkj=self.wpk, m=key)
        self.C = C
        self.W = W

    def Patient_SSHE_encrypt_feature_vector(self):
        encrypted_vector_list = []
        vector_sq_sum = 0
        for e in self.feature_vector:
            encrypted_e = E(m=e, key=self.key)
            encrypted_vector_list.append(encrypted_e)
            vector_sq_sum = gmpy2.add(vector_sq_sum, gmpy2.mul(e, e))
        self.encrypted_feature_vector = encrypted_vector_list
        self.vector_sq_sum = vector_sq_sum

    def Patient_calc_sim(self, sim):
        import global_variable as gl
        self.sim = sim - gmpy2.t_mod_2exp(self.vector_sq_sum, gl.delta1)


    def Patient_Enc(self, wpk_hmc, v):
        return Enc(pkj=self.wpk, m=v), Enc(pkj=wpk_hmc, m=v)


