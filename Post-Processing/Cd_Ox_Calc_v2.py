# filename: Cd_Ox_Calc_v2.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import CoolProp.CoolProp as CP
import os
import glob
from scipy.signal import savgol_filter

from scipy import signal
from Cd_and_dP_avg_calc import Cd_and_dP_avg_calc
from MBV_start_end_idx import MBV_start_end_idx


"""
vvv FILE INPUTS vvv
"""
# Path to test data file
file = r"C:\Users\franc\Downloads\CPL - Post Processing Scripts\All Coldflow and Hotfire Data\TD_CF_1.parquet"

"""
^^^ FILE INPUTS ^^^
"""


"""
vvv ENGINE INPUTS vvv
"""
d_orifice_ox = 0.1285*0.0254
A_ox = 6 * (math.pi/4)*d_orifice_ox**2
d_orifice_ox = 0.040*0.0254        # [m]
A_ox = 0.00000281070638655      # [m2]
A_ox = 21 * (math.pi / 4) * d_orifice_ox**2      # [m2]
A_ox = A_ox + ((math.pi/4) * (0.041*0.0254)**2)
fluid_ox = 'CO2'
# fluid_ox = 'N2O'
A_ox = 0.00000281070638655
p_ox_loss = 105*6894.76

A_ox = 6 * (math.pi/4)*(0.1285*0.0254)**2
A_ox = 0.00001787714657

"""
^^^ ENGINE INPUTS ^^^
"""


"""
vvv WINDOW INPUTS vvv
"""
start_idx = 4220
end_idx = 4289
"""
^^^ WINDOW INPUTS ^^^
"""





print(start_idx)
print(end_idx)
psi_to_pa = 6894.75729
df = pd.read_parquet(file)


df['timestamp'] = pd.to_datetime(df['timestamp'])
df = df.sort_values('timestamp')
df['time_sec'] = (df['timestamp'] - df['timestamp'].iloc[0]).dt.total_seconds()
dt_sample = df['time_sec'].iloc[1]-df['time_sec'].iloc[0]
print(f"dt = {dt_sample} [s]")
print(f"f_s = {1/dt_sample} [Hz]")


# ------------------------------------------------------------
# 1. RAW LOAD CELL (lb → kg)
# ------------------------------------------------------------
print("v")
print("WE STARTING NOW BOYS")
print("v")
# print(f"Max Thrust = {max(df['LC_TOTAL'].values):.2f} [lb]")

t = df['time_sec'].values
m_ox_tank = df['LC_ox_tank'].values * 0.45359237
m_ox_tank_d = m_ox_tank[start_idx:end_idx:5]

t_d = t[start_idx:end_idx:5]

N_Cd = int(end_idx-start_idx)
split = int(N_Cd/5)
t_Cd = t[start_idx:end_idx]
m_ox_tank_Cd = m_ox_tank[start_idx:end_idx]
p1_ox = np.clip(df['PT_ox_tank'].values * psi_to_pa, 0, None)
p2_ox = np.clip(df['PT_chamber'].values * psi_to_pa, 0, None)
p1_ox_Cd = p1_ox[start_idx:end_idx]
p1_ox_Cd = p1_ox_Cd - p_ox_loss
p2_ox_Cd = p2_ox[start_idx:end_idx]
dP_ox_Cd = p1_ox_Cd - p2_ox_Cd
rho_Cd = CP.PropsSI('D', 'P', p1_ox_Cd+101325, 'Q', 0, fluid_ox)
rhodP = rho_Cd * dP_ox_Cd
sqrtrhodP = np.sqrt(rhodP)

split = len(m_ox_tank_d)
Cd_C = np.zeros(split)
dP_C = np.zeros(split)
p1_ox_Cd_split = np.array_split(p1_ox_Cd, split)
p2_ox_Cd_split = np.array_split(p2_ox_Cd, split)
rho_Cd_split = np.array_split(rho_Cd, split)
t_Cd_split = np.array_split(t_Cd,split)
for i in range(split-1):
    Cd_C[i], dP_C[i] = Cd_and_dP_avg_calc(t_Cd_split[i], p1_ox_Cd_split[i], p2_ox_Cd_split[i], m_ox_tank_d[i]-m_ox_tank_d[i+1], rho_Cd[i], A_ox)



plt.figure()
plt.scatter(dP_C*(1/6894.76), Cd_C, label='Raw LC_Ox data')
plt.title('Ox Cd vs. dP')
plt.ylabel('Cd')
plt.xlabel('dP [psi]')
plt.legend()

plt.figure()
plt.scatter(t_d, Cd_C, label='Raw LC_Ox Data')
plt.title('Ox Cd vs. t')
plt.ylabel('Cd')
plt.xlabel('t [s]')
plt.legend()
plt.show()

area_dP = np.trapz(Cd_C, x=dP_C)
area_t = np.trapz(Cd_C, x=t_d)
Cd_avg_dP = (1/(dP_C[-1]-dP_C[0]))*area_dP
Cd_avg_t = (1/(t_d[-1]-t_d[0]))*area_t

print(f"Cd_avg_dP = {Cd_avg_dP}")
print(f"Cd_avg_t = {Cd_avg_t}")


# ============================================================
# SHOW
# ============================================================
plt.show()
