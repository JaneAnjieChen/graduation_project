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


doctor_count = 1
feature_vector = hmc.feature_string_to_vector(feature_string=" ")
# doctors.append(Doctor(index=i, feature_string="", feature_vector=feature_vector))
encrypted_feature_vector, vector_sq_sum = hmc.HMC_SSHE_encrypt_feature_vector(feature_vector)
# hmc send to drp
drp.storage(info_index=0, info=(encrypted_feature_vector, vector_sq_sum),
                IsDoctorInfo=True)


patient_count = 200
patients = []
# 公开的
patients_wpk = {}
for i in range(patient_count):
    tmp = Patient(index=i)
    patients.append(tmp)
    patients_wpk[patient_count] = tmp.publish()
    (v1_C, v1_W), (v_hmc_C, v_hmc_W) = tmp.Patient_Enc(wpk_hmc=wpk_hmc, v=4)
    ehr.storage(doctor_index=0, patient_index=i, v=((v1_C, v1_W), (v_hmc_C, v_hmc_W)))

ehr.extract_v_from_record()
csp_indices = [1, 2, 3, 4, 5, 6, 7]

# 医生评估开始
start = time.clock()

product, count = ehr.mul_v(doctor_index=0)
r_product = drp.calc_r_product_to_csp(product=product, count=count, wpk_hmc=wpk_hmc)

info = [(r_product, 0)]

# drp send to CSPs
for csp in CSPs:
    if csp.index in csp_indices:
        csp.storage(info=info)
        # PDec
        csp.partial_decryption()
        # send to the chosen CSP
        chosen_CSP.storage(info={csp.index: csp.shares}, IsShares=True)
chosen_CSP.combine_shares(csp_indices)
# print(chosen_CSP.result)

doctor_value = drp.calc_doctor_value(m=chosen_CSP.result)

end = time.clock()

print('{}个评分，评估医生，用时{}'.format(patient_count, end-start))

print(doctor_value)
