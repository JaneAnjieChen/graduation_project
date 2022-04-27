from Classes import HMC, CSP, chosen_CSP, DRP, Doctor, Patient, EHR
import time
import numpy as np
from self_defined import dot
import gmpy2

# 1. 系统初始化阶段
print('-'*50 + 'System initializing'+'-'*50)
# start = time.clock()
hmc = HMC()

# 公开的
wpk_hmc, (C_k_hmc, W_k_hmc) = hmc.publish()
# print(hmc.publish())
# pk = wpk_hmc[:2]

import global_variable as gl

hmc.HMC_SkeySplit()
hmc.HMC_calc_B_comb_inverse()
CSPs = []
chosen_CSP = chosen_CSP(index=gl.n, B_comb_inverts=hmc.B_comb_inverts)
for i in range(gl.n):
    if i+1 < gl.n:
        CSPs.append(CSP(index=i, psk=hmc.psk_list[i]))

drp = DRP()
# 公开的
wpk_drp = drp.publish()
k = gl.k

# end = time.clock()
# print(end-start)
ehr = EHR()

print('-'*50 + 'System initialized'+'-'*50)

# 2. 用户注册阶段
print('='*50 + 'Doctor Register Phase' + '='*50)
doctors = []
doctor_count = 0
# d = []
flag = True
while flag:
    flag_str = input('More doctors to be registered (y/n)?')
    if flag_str == 'n':
        flag = False
        break
    elif flag_str == 'y':
        # feature_string 读入
        feature_string = input('Hi doctor, your introduction please (name, hospital, specialty, rank)?')

        # 只有验证通过的医生，在系统中才可以有编号和特征向量
        if hmc.Doctor_verification(feature_string=feature_string):
            doctor_count += 1

            print('Hi doctor, your are now verified.')

            feature_vector = hmc.feature_string_to_vector(feature_string=feature_string)
            doctors.append(Doctor(index=doctor_count-1, feature_string=feature_string,
                                  feature_vector=feature_vector))

            print('Hi doctor, your index is {}'.format(doctor_count-1))

            # start = time.clock()
            encrypted_feature_vector, vector_sq_sum = hmc.HMC_SSHE_encrypt_feature_vector(feature_vector)
            # end = time.clock()
            # print('医生注册用时：{}'.format(end-start))
            # d.append(end-start)

            # hmc send to drp
            drp.storage(info_index=doctor_count-1, info=(encrypted_feature_vector, vector_sq_sum),
                        IsDoctorInfo=True)

            print('Hi doctor {}, you are now registered.'.format(doctor_count-1))

# print(np.average(d))


print('='*50 + 'Patient Register Phase' + '='*50)
patients = []
flag = True
patient_count = 0
# 公开的
patients_wpk = {}

# p = []
while flag:
    patient_count += 1
    flag_str = input('More patients to be registered (y/n)?')
    if flag_str == 'n':
        flag = False
        break
    elif flag_str == 'y':

        # start = time.clock()

        tmp = Patient(index=patient_count-1)

        # end = time.clock()
        # print('患者注册用时：{}'.format(end-start))
        # p.append(end-start)

        patients.append(tmp)
        patients_wpk[patient_count] = tmp.publish()
        print('Hi patient, you are now registered. Your index is {}'.format(patient_count-1))

# print(np.average(p))

# 3. 医生推荐阶段
# s = []
flag = True
while flag:
    flag_str = input('Hello patient, search (y/n)?')
    if flag_str == 'n':
        flag = False
        break
    elif flag_str == 'y':
        print('~'*50 + 'Query begins' + '~'*50)
        patient_index = int(input('Hello patient, your index please?'))

        # patient 数据库里找到index
        for patient in patients:
            if patient.index == patient_index:
                patient_is_requesting = patient
                break

        feature_string = input('Hello patient {}, your request is? '.format(patient_index))
        patient_sim = int(input('Hello patient, set your sim threhold please?(Advised: 0 or 1)'))
        patient_is_requesting.feature_string = feature_string
        patient_is_requesting.feature_string_to_vector()
        # start = time.clock()
        patient_is_requesting.Patient_generate_SSHE_key()
        patient_is_requesting.Patient_SSHE_encrypt_feature_vector()
        patient_is_requesting.Patient_calc_sim(sim=patient_sim)
        # end = time.clock()
        # print(end-start)

        # patient send to drp
        drp.storage(info_index=patient_index, info=(patient_is_requesting.encrypted_feature_vector,
                                                    patient_is_requesting.sim,
                                                    (patient_is_requesting.C, patient_is_requesting.W)),
                    IsPatientInfo=True)

        print('~'*50 + 'Similarity calculation begins' + '~'*50)

        # start = time.clock()

        drp.calc_r_patient_key(patient_index=patient_index)

        info = [(drp.r_C_k_patient, drp.W_k_patient), (C_k_hmc, W_k_hmc)]

        csp_indices = [1, 2, 3, 4, 5, 6, 7]

        # drp send to CSPs
        for csp in CSPs:
            if csp.index in csp_indices:
                csp.storage(info=info)

        for csp in CSPs:
            if csp.index in csp_indices:
                # PDec
                csp.partial_decryption()

                # send to the chosen CSP
                chosen_CSP.storage(info={csp.index: csp.shares}, IsShares=True)

        # hmc.HMC_calc_B_comb_inverse(csp_indices)
        # # send to chosen csp
        # chosen_CSP.storage(info=hmc.B_comb_inverse, IsB_comb_inverse=True)

        # Comb
        chosen_CSP.combine_shares(csp_indices)

        # # check: passed
        # if chosen_CSP.result[0] == gmpy2.mul(patient_is_requesting.key, drp.r):
        #     print('result 0 correct.')
        # else:
        #     print('result 0 incorrect.')
        #
        # if chosen_CSP.result[1] == hmc.key:
        #     print('result 1 correct.')
        # else:
        #     print('result 1 incorrect.')
        #     print(chosen_CSP.result[1])
        #     print(hmc.key)


        # 计算
        chosen_CSP.calc_result_to_drp()

        # 发送给drp
        drp.storage(info_index=None, info=chosen_CSP.k_dec1, k_dec1=True)

        # drp计算相似度 fi
        drp.calc_sim_values_and_recommend_doctors(patient_index=patient_index)

        # # check: passed
        # import global_variable as gl
        # for doctor_index, (doctor_vector, doctor_sq_sum) in drp.doctor_info.items():
        #     print(patient_is_requesting.feature_vector)
        #     # 只有一个医生
        #     print(doctors[0].feature_vector)
        #     dot_product = dot(patient_is_requesting.feature_vector, doctors[0].feature_vector)
        #     dot_product_mod_psai = gmpy2.t_mod(dot_product, gl.psai)
        #     print(dot_product_mod_psai)
        #     print(gl.delta1)
        #     print('真实的明文点乘：', gmpy2.t_mod_2exp(dot_product_mod_psai, gl.delta1))
        # end = time.clock()
        # s.append(end-start)
        # print('相似度计算用时：{}'.format(end-start))
        print('For patient {}, the recommended doctor(s) is/are {}'.format(patient_index, drp.recommend_doctors))

        chosen_doctor_index = int(input('Hi patient {}, which doctor will you choose?'.format(patient_index)))

        if chosen_doctor_index in drp.recommend_doctors:
            print('='*50 + "医生{}向患者{}提供医疗服务".format(chosen_doctor_index, patient_index) + '='*50)
            v = int(input("你好，患者{}，请提交本次医疗服务的反馈评分（0|1|2|3|4|5）：".format(patient_index)))
            # start = time.clock()
            (v1, v_hmc) = patient_is_requesting.Patient_Enc(wpk_hmc=wpk_hmc, v=v)
            # end = time.clock()
            # print(end-start)
            ehr.storage(doctor_index=chosen_doctor_index, patient_index=patient_index, v=(v1, v_hmc))
# print(np.average(s))

# 4. 医生评估阶段
#

