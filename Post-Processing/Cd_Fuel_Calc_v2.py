# filename: Cd_Fuel_Calc_v2.py

import numpy as np
import pandas as pd
import math
import os
import CoolProp.CoolProp as CP
from Cd_and_dP_avg_calc import Cd_and_dP_avg_calc


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
# Injector Area for Fuel [m^2]
d_orifice_inj = 0.021*0.0254
d_orifice_film = 0.021*0.0254
A_film = 7 * (math.pi / 4) * d_orifice_film**2
A_inj = 17 * (math.pi / 4) * d_orifice_inj**2
A_fuel = A_film+A_inj

# Mass Discharged [kg]
m_fuel_discharged = 3.930*0.45359237

# Fluid Density [m^3]
rho_fuel = 1000

# Feedline loss from tank to injector [psi]
P_loss_fuel = 15
"""
^^^ ENGINE INPUTS ^^^
"""


"""
vvv WINDOW INPUTS vvv
"""
start_idx = 4623
end_idx = 4743
"""
^^^ WINDOW INPUTS ^^^
"""

print(start_idx)
print(end_idx)

psi_to_pa = 6894.75729

# ============================================================
# LOAD DATA
# ============================================================
df = pd.read_parquet(file)

df['timestamp'] = pd.to_datetime(df['timestamp'])
df = df.sort_values('timestamp')
df['time_sec'] = (df['timestamp'] - df['timestamp'].iloc[0]).dt.total_seconds()


# ============================================================
# Fuel CD
# ============================================================

print("v")
print("WE STARTING NOW BOYS")
print("v")
print(f"Max Thrust = {max(df['LC_TOTAL'].values):.2f} [lb]")

t = df['time_sec'].values

t_Cd = t[start_idx:end_idx]
p1_fuel = (df['PT_fuel_tank'].values * psi_to_pa)
p2_fuel = (df['PT_chamber'].values * psi_to_pa)
p1_fuel_Cd = p1_fuel[start_idx:end_idx]
p1_fuel_Cd = p1_fuel[start_idx:end_idx] - (P_loss_fuel*psi_to_pa*np.ones(len(p1_fuel_Cd)))
p2_fuel_Cd = p2_fuel[start_idx:end_idx]

rho_Cd = CP.PropsSI('D', 'P', p1_fuel_Cd - (P_loss_fuel*psi_to_pa*np.ones(len(p1_fuel_Cd))), 'T', 298, 'H2O')

Cd_fuel, dP_fuel = Cd_and_dP_avg_calc(t_Cd, p1_fuel_Cd, p2_fuel_Cd, m_fuel_discharged, rho_Cd, A_fuel)
print(f"dP_fuel = {dP_fuel/6894.76:.2f} [psi]")
print(f"Cd_fuel = {Cd_fuel}")
