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

# doctors = []
doctor_count = 200
for i in range(doctor_count):
    feature_vector = hmc.feature_string_to_vector(feature_string=" ")
    # doctors.append(Doctor(index=i, feature_string="", feature_vector=feature_vector))
    encrypted_feature_vector, vector_sq_sum = hmc.HMC_SSHE_encrypt_feature_vector(feature_vector)
    # hmc send to drp
    drp.storage(info_index=i, info=(encrypted_feature_vector, vector_sq_sum),
                IsDoctorInfo=True)

print(drp.doctor_info.keys())

patient_count = 1
patients_wpk = {}
patient = Patient(index=0)
patients_wpk[0] = patient.publish()


patient_index = 0
patient_sim = 1
patient.feature_string = " "
patient.feature_string_to_vector()
patient.Patient_generate_SSHE_key()
patient.Patient_SSHE_encrypt_feature_vector()
patient.Patient_calc_sim(sim=patient_sim)
# patient send to drp
drp.storage(info_index=patient_index, info=(patient.encrypted_feature_vector,
                                            patient.sim, (patient.C, patient.W)), IsPatientInfo=True)




print('~'*50 + 'Similarity calculation begins' + '~'*50)

start = time.clock()

drp.calc_r_patient_key(patient_index=patient_index)
info = [(drp.r_C_k_patient, drp.W_k_patient), (C_k_hmc, W_k_hmc)]

csp_indices = [1, 2, 3, 4, 5, 6, 7]
# drp send to CSPs
for csp in CSPs:
    if csp.index in csp_indices:
        csp.storage(info=info)
        # PDec
        csp.partial_decryption()
        # send to the chosen CSP
        chosen_CSP.storage(info={csp.index: csp.shares}, IsShares=True)

chosen_CSP.combine_shares(csp_indices)
chosen_CSP.calc_result_to_drp()
# 发送给drp
drp.storage(info_index=None, info=chosen_CSP.k_dec1, k_dec1=True)
# drp计算相似度 fi
drp.calc_sim_values_and_recommend_doctors(patient_index=patient_index)

end = time.clock()

print('{}个医生相似度计算用时：{}'.format(doctor_count, end-start))
print('For patient {}, the recommended doctor(s) is/are {}'.format(patient_index, drp.recommend_doctors))
