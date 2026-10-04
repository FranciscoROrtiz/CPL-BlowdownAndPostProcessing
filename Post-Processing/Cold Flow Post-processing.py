# filename: Cold Flow Post-processing.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import CoolProp.CoolProp as CP
from scipy.signal import savgol_filter

##### Initial Parameters #####
Cd = 0.6        # Discharge Coefficient

d_orifice_fuel = 0.0005334          # Fuel orifice diameter [m]; #75 drill bit
A_orifice_fuel = ( math.pi / 4 ) * (d_orifice_fuel ** 2)     # Fuel orifice area [m^2]
n_fuel = 17     # Number of fuel orifices
A_fuel = n_fuel * A_orifice_fuel        # Total Fuel Orifice Area

T = 298     # Temperature [K]
rho = CP.PropsSI('D', 'P', p1_list[0], 'T', T, fluid)      # Density [kg/m^3]
#delta_t_LC =        # CYCLE time in Revel code


##### PT Mass Flow Rate #####
mdot_PT = []
for i in range(len(p1_list)):
    p1 = p1_list[i] * 6894.75729        # Inlet pressure [Pa]; [psi] * 6894.75729 = [Pa]
    p2 = p2_list[i] * 6894.75729        # Outlet pressure [Pa]; [psi] * 6894.75729 = [Pa]

    if p1 <= p2:
        mdot_PT.append(0)
        continue        # "continue" skips the rest of the iteration and moves to next loop

    mdot = (Cd * A_fuel) * ( (2 * rho * (p1 - p2)) ** 0.5 )  # Calculate mass flow rate from PTs with SPI model
    mdot_PT.append(mdot)


##### Load Cell Mass Flow Rate #####
mdot_LC = []
for i in range(1, len(load_cell_data)):
    mass_previous = load_cell_data[i-1]
    mass_current = load_cell_data[i]
    mdot = (mass_previous - mass_current) / delta_t_LC      # mdot should be greater than 0 (positive)
    mdot_LC.append(mdot)


##### Error #####

mdot_LC = np.array(mdot_LC)
mdot_PT = np.array(mdot_PT)
error = mdot_LC - mdot_PT
error_percent = 100 * (mdot_PT - mdot_LC) / mdot_LC

plt.plot(time[1:], error_percent)
plt.ylabel("Percent Error")
plt.xlabel("Time")
plt.grid()
plt.show()
