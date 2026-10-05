# filename: PP_Ox_EM_Fuel_SPI_v4.py
# created: 04.22.26
# last edited: 04.23.26
# authors: Francisco Ortiz, Diego Ortiz
# purpose: use engine design paramters and propellant conditions to estimate initial propellant loading conditions required and resulatant propellant blowdown time history



import json
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import time
import pandas as pd
import math as math
import cantera as ct
from scipy.optimize import brentq
import ast

from generatePlot import generatePlot
from def_combustion_gas_new import def_combustion_gas_new
from calc_chamber_pressure import calc_chamber_pressure
from area_mach_residual import area_mach_residual
from thrust import thrust
from mdot_SPI import mdot_SPI
from AdiabaticBlowdown import AdiabaticBlowdown
from calculate_orifice_Areas_Fuel import calculate_orifice_Areas_Fuel
from fuel_time_step import fuel_time_step
from define_conditions import define_conditions
from calculate_mass_flux_rate import calculate_mass_flux_rate
from calculate_injector_Area import calculate_injector_Area
from ox_temp_calc import ox_temp_calc
from ox_time_step import ox_time_step
from N2O_ullage import N2O_ullage_TVm
from N2O_ullage import N2O_ullage_PVm
from N2O_tank_mass import N2O_tank_mass_TVml
from N2O_tank_mass import N2O_tank_mass_PVml
from N2O_tank_mass import N2O_liquid_mass_PVm
from save_data import save_data_json_single
"""
Output =
[
t,
Pc,
OF,
phi,
P_pressurant,
P_injector,
mdot_fuel_comb,
mdot_fuel_film,
m_fuel,
V_fuel,
V_pressurant,
T_ox,
P_ox,
mdot_ox,
m_ox,
m_ox_vapor,
m_ox_liquid,
U_tot,
Thrust
]
"""

def blowdown_solve(Input, Files):

    # Assign variables from function inpuits to readable variables
    Cd_ox = Input['Cd_ox']
    A_ox = Input['A_ox']
    P_loss_ox = Input['P_loss_ox']
    P_stiffness_ox = Input['P_stiffness_ox']
    m_ox_0 = Input['m_ox_0']
    Substance = Input['Substance']
    P_ox_0 = Input['P_ox_0']
    rho_l_tab = Files['rho_l_tab']
    rho_v_tab = Files['rho_v_tab']
    u_l_tab = Files['u_l_tab']
    u_v_tab = Files['u_v_tab']
    T_grid = Files['T_grid']
    V_tank_ox = Input['V_tank_ox']

    Cd_fuel = Input['Cd_fuel']
    A_inj = Input['A_inj']
    A_film = Input['A_film']
    P_stiffness_fuel = Input['P_stiffness_fuel']
    P_loss_fuel = Input['P_loss_fuel']
    rho_fuel = Input['rho_fuel']
    gamma = Input['gamma']
    P_pressurant_0 = Input['P_pressurant_0']
    V_tank_fuel = Input['V_tank_fuel']
    Ullage_fuel = Input['Ullage_fuel']

    t_burn = Input['t_burn']
    Pc_0 = Input['Pc_0']
    A_th = Input['A_th']
    A_e = Input['A_e']
    OF_st = Input['OF_st']
    n_cstar = Input['n_cstar']
    phi_tab = Files['phi_tab']
    Tad_tab = Files['Tad_tab']
    gamma_tab = Files['gamma_tab']
    R_tab = Files['R_tab']


    # Set timestep
    dt = 0.01      # [s]

    # Calculate number of points in time-history vectors
    N = int((t_burn*10 / dt) + 1)

    # Initizalize arrays
    t = np.zeros(N)
    Pc = np.zeros(N)

    T_ox = np.zeros(N)
    P_ox = np.zeros(N)
    m_ox = np.zeros(N)
    m_ox_vapor = np.zeros(N)
    m_ox_liquid = np.zeros(N)
    mdot_ox = np.zeros(N)
    U_tot = np.zeros(N)

    P_pressurant = np.zeros(N)
    P_injector = np.zeros(N)
    m_fuel = np.zeros(N)
    mdot_fuel_comb = np.zeros(N)
    mdot_fuel_film = np.zeros(N)
    V_fuel = np.zeros(N)
    V_pressurant = np.zeros(N)

    OF = np.zeros(N)
    phi = np.zeros(N)
    Thrust = np.zeros(N)

    # Initial Liquid and Vapor Properties
    T_ox_0 = CP.PropsSI('T', 'P', P_ox_0, 'Q', 1, Substance)    # Temperature of the N2O [K]
    m_ox_liquid_0 = N2O_liquid_mass_PVm(P_ox_0, V_tank_ox, m_ox_0)
    P0 = CP.PropsSI('P', 'T', T_ox_0, 'Q', 1, Substance)    # Pressure of the N2O [Pa]
    rho_l = CP.PropsSI('D', 'T', T_ox_0, 'Q', 0, Substance)  # Density of liquid
    rho_v = CP.PropsSI('D', 'T', T_ox_0, 'Q', 1, Substance)  # Density of vapor
    u_l = CP.PropsSI('U', 'T', T_ox_0, 'Q', 0, Substance)  # Internal energy of liquid
    u_v = CP.PropsSI('U', 'T', T_ox_0, 'Q', 1, Substance)  # Internal energy of vapor
    # Initial Volumes and Masses
    V_l = m_ox_liquid_0 / rho_l
    V_v = V_tank_ox - V_l
    # V_v = V_l * (Ullage_ox/(1-Ullage_ox))
    m_v = rho_v * V_v

    # V_tank_ox = V_l + V_v
    print(f"V_tank_ox = {V_tank_ox*1000:.2f} [L]")

    mtot0 = m_v + m_ox_liquid_0
    x = m_v / mtot0

    Pc[0] = Pc_0

    # Initial Total Energy of Oxidizer
    U_tot[0] = (m_v * u_v) + (m_ox_liquid_0 * u_l)

    # Set Initial Values
    T_ox[0] = T_ox_0
    P_ox[0] = P_ox_0
    m_ox[0] = mtot0
    m_ox_vapor[0] = m_v
    m_ox_liquid[0] = m_ox_liquid_0

    P_pressurant[0] = P_pressurant_0
    P_injector[0] = P_pressurant_0 - P_loss_fuel
    V_pressurant[0] = Ullage_fuel * V_tank_fuel
    V_fuel[0] = V_tank_fuel - V_pressurant[0]
    m_fuel[0] = rho_fuel * V_fuel[0]

    Stiffness_warning = 0
    i = 0
    while True:
        if i*dt == 0:
            print("0 seconds")

        # print(f"Iteration: {i}")
        t[i + 1] = t[i] + dt
        T_ox[i + 1], P_ox[i + 1], mdot_ox[i], m_ox[i + 1], m_ox_vapor[i + 1], m_ox_liquid[i + 1], U_tot[i + 1] = ox_time_step(P_loss_ox, Cd_ox, A_ox, Substance, V_tank_ox, dt, P_ox[i], m_ox[i], T_ox[i], U_tot[i], Pc[i], rho_l_tab, rho_v_tab, u_l_tab, u_v_tab, T_grid)
        P_pressurant[i + 1], P_injector[i + 1], mdot_fuel_comb[i], mdot_fuel_film[i], m_fuel[i + 1], V_fuel[i + 1], V_pressurant[i + 1] = fuel_time_step(P_loss_fuel, Cd_fuel, A_inj, A_film, rho_fuel, dt, P_pressurant[i], m_fuel[i], V_fuel[i], V_pressurant[i], gamma, Pc[i])

        OF[i] = mdot_ox[i] / mdot_fuel_comb[i]
        phi[i] = OF_st / OF[i]

        Pc[i + 1] = calc_chamber_pressure(phi[i], (mdot_ox[i]+mdot_fuel_comb[i]), A_th, phi_tab, Tad_tab, gamma_tab, R_tab, n_cstar)

        Thrust[i + 1] = thrust(Pc[i+1], phi[i], (mdot_ox[i]+mdot_fuel_comb[i]), A_e, A_th, phi_tab, Tad_tab, gamma_tab, R_tab)

        # Ox Thresholds
        if T_ox[i + 1] <= 182.23:  # Stop further computation if cutoff is reached
            print(f"OXIDIZER: Temperature hit cutoff (182.23[K]) at iteration {i}")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            return Output
        # if m_ox_liquid[i + 1] <= 0 or (m_ox_liquid[i+1]>m_ox_liquid[i]):
        if m_ox_liquid[i + 1] <= 0:
            print(f"OXIDIZER: Liquid mass depleted at iteration {i}: {m_ox_liquid[i + 1]:.6f} kg & {m_ox_liquid[i]}")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            return Output
        if P_ox[i+1] < (Pc[i]+(Pc[i]*P_stiffness_ox)+P_loss_ox):
            print(f"OXIDIZER: Pressure below stiffness condition at iteration {i}: {P_ox[i+1]/6894.76:.2f}[psi]")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            return Output

        # Fuel Thresholds
        if m_fuel[i + 1] <= 0:
            print(f"\nERROR:")
            print(f"FUEL: Liquid mass depleted at iteration {i}: {m_fuel[i + 1]:.6f} kg\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            return Output
        # Warning that Pressure is at stiffness condition
        if ((P_pressurant[i+1] - P_loss_fuel - Pc[i] - (Pc[i]*P_stiffness_fuel)) <= 1) and (Stiffness_warning == 0):
            Stiffness_warning = 1
            print(f"\nWARNING:")
            print(f"FUEL: Pressure is at Stiffness condition at iteration {i}: {t[i+1]:.3f} [s]\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            return Output
        # Stop if Pressure of fuel is less than or equal to chamber pressure
        del_P = P_pressurant[i+1] - P_loss_fuel - Pc[i]
        if del_P <= 1:
            print(f"\nERROR:")
            print(f"FUEL: Fuel Pressure went below minimum injector pressure limit at iteration {i}: {P_pressurant[i + 1]*0.000145038:.6f} psi\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            return Output

        i += 1




    dP = P_pressurant[-1] - P_loss_fuel - Pc[-1]
    mdot_fuel_comb[-1] = mdot_SPI(Cd_fuel, A_inj, rho_fuel, dP)
    mdot_fuel_film[-1] = mdot_SPI(Cd_fuel, A_film, rho_fuel, dP)
    N2O = define_conditions((P_ox[-1] - P_loss_ox), Pc[-1])
    G = calculate_mass_flux_rate(N2O)
    mdot_ox[-1] = G*Cd_ox*A_ox
    OF[-1] = mdot_ox[-1] / mdot_fuel_comb[-1]
    phi[-1] = OF_st / OF[-1]

    Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
    return Output

def blowdown_solve_set_tburn(Input, Files):
    # Ox=[Cd_ox, A_ox, P_loss_ox, P_stiffness_ox, Ullage_ox, m_ox_liquid_0, Substance, T_ox_0],
    # Fuel=[Cd_fuel, A_inj, A_film, P_stiffness_fuel, P_loss_fuel, rho_fuel, gamma, P_pressurant_0, V_tank_fuel, Ullage_fuel],
    # Engine=[t_burn, Pc_0, A_th, mech, OF_st]

    # Assign variables from function inpuits to readable variables
    Cd_ox = Input['Cd_ox']
    A_ox = Input['A_ox']
    P_loss_ox = Input['P_loss_ox']
    P_stiffness_ox = Input['P_stiffness_ox']
    m_ox_0 = Input['m_ox_0']
    Substance = Input['Substance']
    P_ox_0 = Input['P_ox_0']
    rho_l_tab = Files['rho_l_tab']
    rho_v_tab = Files['rho_v_tab']
    u_l_tab = Files['u_l_tab']
    u_v_tab = Files['u_v_tab']
    T_grid = Files['T_grid']
    V_tank_ox = Input['V_tank_ox']

    Cd_fuel = Input['Cd_fuel']
    A_inj = Input['A_inj']
    A_film = Input['A_film']
    P_stiffness_fuel = Input['P_stiffness_fuel']
    P_loss_fuel = Input['P_loss_fuel']
    rho_fuel = Input['rho_fuel']
    gamma = Input['gamma']
    P_pressurant_0 = Input['P_pressurant_0']
    V_tank_fuel = Input['V_tank_fuel']
    Ullage_fuel = Input['Ullage_fuel']

    t_burn = Input['t_burn']
    Pc_0 = Input['Pc_0']
    A_th = Input['A_th']
    A_e = Input['A_e']
    OF_st = Input['OF_st']
    n_cstar = Input['n_cstar']
    phi_tab = Files['phi_tab']
    Tad_tab = Files['Tad_tab']
    gamma_tab = Files['gamma_tab']
    R_tab = Files['R_tab']

    # Set timestep
    dt = 0.01      # [s]

    # Calculate number of points in time-history vectors
    N = int((t_burn / dt) + 1)

    # Initizalize arrays
    t = np.zeros(N)
    Pc = np.zeros(N)

    T_ox = np.zeros(N)
    P_ox = np.zeros(N)
    m_ox = np.zeros(N)
    m_ox_vapor = np.zeros(N)
    m_ox_liquid = np.zeros(N)
    mdot_ox = np.zeros(N)
    U_tot = np.zeros(N)

    P_pressurant = np.zeros(N)
    P_injector = np.zeros(N)
    m_fuel = np.zeros(N)
    mdot_fuel_comb = np.zeros(N)
    mdot_fuel_film = np.zeros(N)
    V_fuel = np.zeros(N)
    V_pressurant = np.zeros(N)

    OF = np.zeros(N)
    phi = np.zeros(N)
    Thrust = np.zeros(N)

    # Initial Liquid and Vapor Properties
    T_ox_0 = CP.PropsSI('T', 'P', P_ox_0, 'Q', 1, Substance)    # Temperature of the N2O [K]
    m_ox_liquid_0 = N2O_liquid_mass_PVm(P_ox_0, V_tank_ox, m_ox_0)
    rho_l = CP.PropsSI('D', 'T', T_ox_0, 'Q', 0, Substance)  # Density of liquid
    rho_v = CP.PropsSI('D', 'T', T_ox_0, 'Q', 1, Substance)  # Density of vapor
    u_l = CP.PropsSI('U', 'T', T_ox_0, 'Q', 0, Substance)  # Internal energy of liquid
    u_v = CP.PropsSI('U', 'T', T_ox_0, 'Q', 1, Substance)  # Internal energy of vapor
    # Initial Volumes and Masses
    V_l = m_ox_liquid_0 / rho_l
    V_v = V_tank_ox - V_l
    # V_v = V_l * (Ullage_ox/(1-Ullage_ox))
    m_v = rho_v * V_v

    # V_tank_ox = V_l + V_v
    print(f"V_tank_ox = {V_tank_ox*1000:.2f} [L]")

    mtot0 = m_v + m_ox_liquid_0
    x = m_v / mtot0

    Pc[0] = Pc_0

    # Initial Total Energy of Oxidizer
    U_tot[0] = (m_v * u_v) + (m_ox_liquid_0 * u_l)

    # Set Initial Values
    T_ox[0] = T_ox_0
    P_ox[0] = P_ox_0
    m_ox[0] = mtot0
    m_ox_vapor[0] = m_v
    m_ox_liquid[0] = m_ox_liquid_0

    P_pressurant[0] = P_pressurant_0
    P_injector[0] = P_pressurant_0 - P_loss_fuel
    V_pressurant[0] = Ullage_fuel * V_tank_fuel
    V_fuel[0] = V_tank_fuel - V_pressurant[0]
    m_fuel[0] = rho_fuel * V_fuel[0]

    Stiffness_warning = 0
    i = 0
    while (i*dt) < (t_burn):
        if i*dt == 0:
            print("0 seconds")

        # print(f"Iteration: {i}")
        t[i + 1] = t[i] + dt
        T_ox[i + 1], P_ox[i + 1], mdot_ox[i], m_ox[i + 1], m_ox_vapor[i + 1], m_ox_liquid[i + 1], U_tot[i + 1] = ox_time_step(P_loss_ox, Cd_ox, A_ox, Substance, V_tank_ox, dt, P_ox[i], m_ox[i], T_ox[i], U_tot[i], Pc[i], rho_l_tab, rho_v_tab, u_l_tab, u_v_tab, T_grid)
        P_pressurant[i + 1], P_injector[i + 1], mdot_fuel_comb[i], mdot_fuel_film[i], m_fuel[i + 1], V_fuel[i + 1], V_pressurant[i + 1] = fuel_time_step(P_loss_fuel, Cd_fuel, A_inj, A_film, rho_fuel, dt, P_pressurant[i], m_fuel[i], V_fuel[i], V_pressurant[i], gamma, Pc[i])

        OF[i] = mdot_ox[i] / mdot_fuel_comb[i]
        phi[i] = OF_st / OF[i]

        Pc[i + 1] = calc_chamber_pressure(phi[i], (mdot_ox[i]+mdot_fuel_comb[i]), A_th, phi_tab, Tad_tab, gamma_tab, R_tab, n_cstar)

        Thrust[i + 1] = thrust(Pc[i+1], phi[i], (mdot_ox[i]+mdot_fuel_comb[i]), A_e, A_th, phi_tab, Tad_tab, gamma_tab, R_tab)

        # Ox Thresholds
        if T_ox[i + 1] <= 182.23:  # Stop further computation if cutoff is reached
            print(f"OXIDIZER: Temperature hit cutoff (182.23[K]) at iteration {i}")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("ox")
            return Output
        # if m_ox_liquid[i + 1] <= 0 or (m_ox_liquid[i+1]>m_ox_liquid[i]):
        if m_ox_liquid[i + 1] <= 0:
            print(f"OXIDIZER: Liquid mass depleted at iteration {i}: {m_ox_liquid[i + 1]:.6f} kg & {m_ox_liquid[i]}")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("ox")
            return Output
        if P_ox[i+1] < (Pc[i]+(Pc[i]*P_stiffness_ox)+P_loss_ox):
            print(f"OXIDIZER: Pressure below stiffness condition at iteration {i}: {P_ox[i+1]/6894.76:.2f}[psi]")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("ox")
            return Output

        # Fuel Thresholds
        if m_fuel[i + 1] <= 0:
            print(f"\nERROR:")
            print(f"FUEL: Liquid mass depleted at iteration {i}: {m_fuel[i + 1]:.6f} kg\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("fuel")
            return Output
        # Warning that Pressure is at stiffness condition
        if ((P_pressurant[i+1] - P_loss_fuel - Pc[i] - (Pc[i]*P_stiffness_fuel)) <= 1) and (Stiffness_warning == 0) and (i>5):
            Stiffness_warning = 1
            print(f"\nWARNING:")
            print(f"FUEL: Pressure is at Stiffness condition at iteration {i}: {t[i+1]:.3f} [s]\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("fuel")
            return Output
        # Stop if Pressure of fuel is less than or equal to chamber pressure
        del_P = P_pressurant[i+1] - P_loss_fuel - Pc[i]
        if del_P <= 1:
            print(f"\nERROR:")
            print(f"FUEL: Fuel Pressure went below minimum injector pressure limit at iteration {i}: {P_pressurant[i + 1]*0.000145038:.6f} psi\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("fuel")
            return Output

        i += 1



    dP = P_pressurant[-1] - P_loss_fuel - Pc[-1]
    mdot_fuel_comb[-1] = mdot_SPI(Cd_fuel, A_inj, rho_fuel, dP)
    mdot_fuel_film[-1] = mdot_SPI(Cd_fuel, A_film, rho_fuel, dP)
    N2O = define_conditions((P_ox[-1] - P_loss_ox), Pc[-1])
    G = calculate_mass_flux_rate(N2O)
    mdot_ox[-1] = G*Cd_ox*A_ox
    OF[-1] = mdot_ox[-1] / mdot_fuel_comb[-1]
    phi[-1] = OF_st / OF[-1]

    Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
    return Output

def blowdown_solve_tank_loading(Input, Files):
    while True:
        data = blowdown_solve_set_tburn(Input, Files)
        if data[0][-1]>=(Input['t_burn']-0.01):
            break
        elif data[-1] == "ox":
            Input['m_ox_0'] += 0.01
        else:
            Input['Ullage_fuel'] -= 0.01
    return data

def coldflow_solve_set_tburn(Input, Files):
    # Ox=[Cd_ox, A_ox, P_loss_ox, P_stiffness_ox, Ullage_ox, m_ox_liquid_0, Substance, T_ox_0],
    # Fuel=[Cd_fuel, A_inj, A_film, P_stiffness_fuel, P_loss_fuel, rho_fuel, gamma, P_pressurant_0, V_tank_fuel, Ullage_fuel],
    # Engine=[t_burn, Pc_0, A_th, mech, OF_st]

    # Assign variables from function inpuits to readable variables
    Cd_ox = Input['Cd_ox']
    A_ox = Input['A_ox']
    P_loss_ox = Input['P_loss_ox']
    P_stiffness_ox = Input['P_stiffness_ox']
    m_ox_0 = Input['m_ox_0']
    Substance = Input['Substance']
    P_ox_0 = Input['P_ox_0']
    rho_l_tab = Files['rho_l_tab']
    rho_v_tab = Files['rho_v_tab']
    u_l_tab = Files['u_l_tab']
    u_v_tab = Files['u_v_tab']
    T_grid = Files['T_grid']
    V_tank_ox = Input['V_tank_ox']

    Cd_fuel = Input['Cd_fuel']
    A_inj = Input['A_inj']
    A_film = Input['A_film']
    P_stiffness_fuel = Input['P_stiffness_fuel']
    P_loss_fuel = Input['P_loss_fuel']
    rho_fuel = Input['rho_fuel']
    gamma = Input['gamma']
    P_pressurant_0 = Input['P_pressurant_0']
    V_tank_fuel = Input['V_tank_fuel']
    Ullage_fuel = Input['Ullage_fuel']

    t_burn = Input['t_burn']
    Pc_0 = Input['Pc_0']
    A_th = Input['A_th']
    A_e = Input['A_e']
    OF_st = Input['OF_st']
    n_cstar = Input['n_cstar']
    phi_tab = Files['phi_tab']
    Tad_tab = Files['Tad_tab']
    gamma_tab = Files['gamma_tab']
    R_tab = Files['R_tab']

    # Set timestep
    dt = 0.01      # [s]

    # Calculate number of points in time-history vectors
    N = int((t_burn / dt) + 1)

    # Initizalize arrays
    t = np.zeros(N)
    Pc = np.zeros(N)

    T_ox = np.zeros(N)
    P_ox = np.zeros(N)
    m_ox = np.zeros(N)
    m_ox_vapor = np.zeros(N)
    m_ox_liquid = np.zeros(N)
    mdot_ox = np.zeros(N)
    U_tot = np.zeros(N)

    P_pressurant = np.zeros(N)
    P_injector = np.zeros(N)
    m_fuel = np.zeros(N)
    mdot_fuel_comb = np.zeros(N)
    mdot_fuel_film = np.zeros(N)
    V_fuel = np.zeros(N)
    V_pressurant = np.zeros(N)

    OF = np.zeros(N)
    phi = np.zeros(N)
    Thrust = np.zeros(N)

    # Initial Liquid and Vapor Properties
    T_ox_0 = CP.PropsSI('T', 'P', P_ox_0, 'Q', 1, Substance)    # Temperature of the N2O [K]
    m_ox_liquid_0 = N2O_liquid_mass_PVm(P_ox_0, V_tank_ox, m_ox_0)
    rho_l = CP.PropsSI('D', 'T', T_ox_0, 'Q', 0, Substance)  # Density of liquid
    rho_v = CP.PropsSI('D', 'T', T_ox_0, 'Q', 1, Substance)  # Density of vapor
    u_l = CP.PropsSI('U', 'T', T_ox_0, 'Q', 0, Substance)  # Internal energy of liquid
    u_v = CP.PropsSI('U', 'T', T_ox_0, 'Q', 1, Substance)  # Internal energy of vapor
    # Initial Volumes and Masses
    V_l = m_ox_liquid_0 / rho_l
    V_v = V_tank_ox - V_l
    # V_v = V_l * (Ullage_ox/(1-Ullage_ox))
    m_v = rho_v * V_v

    # V_tank_ox = V_l + V_v
    print(f"V_tank_ox = {V_tank_ox*1000:.2f} [L]")

    mtot0 = m_v + m_ox_liquid_0
    x = m_v / mtot0

    Pc[0] = Pc_0

    # Initial Total Energy of Oxidizer
    U_tot[0] = (m_v * u_v) + (m_ox_liquid_0 * u_l)

    # Set Initial Values
    T_ox[0] = T_ox_0
    P_ox[0] = P_ox_0
    m_ox[0] = mtot0
    m_ox_vapor[0] = m_v
    m_ox_liquid[0] = m_ox_liquid_0

    P_pressurant[0] = P_pressurant_0
    P_injector[0] = P_pressurant_0 - P_loss_fuel
    V_pressurant[0] = Ullage_fuel * V_tank_fuel
    V_fuel[0] = V_tank_fuel - V_pressurant[0]
    m_fuel[0] = rho_fuel * V_fuel[0]

    Stiffness_warning = 0
    i = 0
    while (i*dt) < (t_burn):
        if i*dt == 0:
            print("0 seconds")

        # print(f"Iteration: {i}")
        t[i + 1] = t[i] + dt
        T_ox[i + 1], P_ox[i + 1], mdot_ox[i], m_ox[i + 1], m_ox_vapor[i + 1], m_ox_liquid[i + 1], U_tot[i + 1] = ox_time_step(P_loss_ox, Cd_ox, A_ox, Substance, V_tank_ox, dt, P_ox[i], m_ox[i], T_ox[i], U_tot[i], Pc[i], rho_l_tab, rho_v_tab, u_l_tab, u_v_tab, T_grid)
        P_pressurant[i + 1], P_injector[i + 1], mdot_fuel_comb[i], mdot_fuel_film[i], m_fuel[i + 1], V_fuel[i + 1], V_pressurant[i + 1] = fuel_time_step(P_loss_fuel, Cd_fuel, A_inj, A_film, rho_fuel, dt, P_pressurant[i], m_fuel[i], V_fuel[i], V_pressurant[i], gamma, Pc[i])

        OF[i] = mdot_ox[i] / mdot_fuel_comb[i]
        phi[i] = OF_st / OF[i]

        Pc[i + 1] = 101325

        Thrust[i + 1] = thrust(Pc[i+1], phi[i], (mdot_ox[i]+mdot_fuel_comb[i]), A_e, A_th, phi_tab, Tad_tab, gamma_tab, R_tab)

        # Ox Thresholds
        if T_ox[i + 1] <= 182.23:  # Stop further computation if cutoff is reached
            print(f"OXIDIZER: Temperature hit cutoff (182.23[K]) at iteration {i}")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("ox")
            return Output
        # if m_ox_liquid[i + 1] <= 0 or (m_ox_liquid[i+1]>m_ox_liquid[i]):
        if m_ox_liquid[i + 1] <= 0:
            print(f"OXIDIZER: Liquid mass depleted at iteration {i}: {m_ox_liquid[i + 1]:.6f} kg & {m_ox_liquid[i]}")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("ox")
            return Output
        if P_ox[i+1] < (Pc[i]+(Pc[i]*P_stiffness_ox)+P_loss_ox):
            print(f"OXIDIZER: Pressure below stiffness condition at iteration {i}: {P_ox[i+1]/6894.76:.2f}[psi]")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("ox")
            return Output

        # Fuel Thresholds
        if m_fuel[i + 1] <= 0:
            print(f"\nERROR:")
            print(f"FUEL: Liquid mass depleted at iteration {i}: {m_fuel[i + 1]:.6f} kg\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("fuel")
            return Output
        # Warning that Pressure is at stiffness condition
        if ((P_pressurant[i+1] - P_loss_fuel - Pc[i] - (Pc[i]*P_stiffness_fuel)) <= 1) and (Stiffness_warning == 0):
            Stiffness_warning = 1
            print(f"\nWARNING:")
            print(f"FUEL: Pressure is at Stiffness condition at iteration {i}: {t[i+1]:.3f} [s]\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("fuel")
            return Output
        # Stop if Pressure of fuel is less than or equal to chamber pressure
        del_P = P_pressurant[i+1] - P_loss_fuel - Pc[i]
        if del_P <= 1:
            print(f"\nERROR:")
            print(f"FUEL: Fuel Pressure went below minimum injector pressure limit at iteration {i}: {P_pressurant[i + 1]*0.000145038:.6f} psi\n")
            Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
            Output = [item[:i+1] for item in Output]
            Output.append("fuel")
            return Output

        i += 1



    dP = P_pressurant[-1] - P_loss_fuel - Pc[-1]
    mdot_fuel_comb[-1] = mdot_SPI(Cd_fuel, A_inj, rho_fuel, dP)
    mdot_fuel_film[-1] = mdot_SPI(Cd_fuel, A_film, rho_fuel, dP)
    N2O = define_conditions((P_ox[-1] - P_loss_ox), Pc[-1])
    G = calculate_mass_flux_rate(N2O)
    mdot_ox[-1] = G*Cd_ox*A_ox
    OF[-1] = mdot_ox[-1] / mdot_fuel_comb[-1]
    phi[-1] = OF_st / OF[-1]

    Output = [t, Pc, OF, phi, P_pressurant, P_injector, mdot_fuel_comb, mdot_fuel_film, m_fuel, V_fuel, V_pressurant, T_ox, P_ox, mdot_ox, m_ox, m_ox_vapor, m_ox_liquid, U_tot, Thrust]
    return Output
