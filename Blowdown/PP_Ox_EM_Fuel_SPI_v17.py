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
from save_data import save_data_json_single
from save_data import save_data_json
from blowdown_MK2 import blowdown_solve
from blowdown_MK2 import blowdown_solve_set_tburn
from blowdown_MK2 import blowdown_solve_tank_loading

def main():

    """
    vvv INPUTS vvv
    """
    Input_file_path = r'EngineInput_TD_HF_1.txt'

    rho_l_tab_json_file_path = "rho_l_CoolPropProperties_N2O_dT0.001_T1182.23_T2305.json"
    rho_v_tab_json_file_path = "rho_v_CoolPropProperties_N2O_dT0.001_T1182.23_T2305.json"
    u_l_tab_json_file_path = "u_l_CoolPropProperties_N2O_dT0.001_T1182.23_T2305.json"
    u_v_tab_json_file_path = "u_v_CoolPropProperties_N2O_dT0.001_T1182.23_T2305.json"
    Tad_tab_json_file_path = "Tad_IPA_N2O_EquilibrateProperties_dphi0.010_phi01.000.json"
    gamma_tab_json_file_path = "gamma_IPA_N2O_EquilibrateProperties_dphi0.010_phi01.000.json"
    R_tab_json_file_path = "R_IPA_N2O_EquilibrateProperties_dphi0.010_phi01.000.json"
    T_grid_json_file_path = "T_grid.json"
    phi_tab_json_file_path = "phi_tab.json"
    """
    ^^^ INPUTS ^^^
    """



    """
    vvv LOAD RELEVANT DATA (DON'T DELETE)vvv
    """
    with open(Input_file_path, 'r') as file:
        data = file.read()
        Input = ast.literal_eval(data)

    with open(rho_l_tab_json_file_path, 'r') as file:
        rho_l_tab = json.load(file)


    with open(rho_v_tab_json_file_path, 'r') as file:
        rho_v_tab = json.load(file)


    with open(u_l_tab_json_file_path, 'r') as file:
        u_l_tab = json.load(file)


    with open(u_v_tab_json_file_path, 'r') as file:
        u_v_tab = json.load(file)

    with open(Tad_tab_json_file_path, 'r') as file:
        Tad_tab = json.load(file)

    with open(gamma_tab_json_file_path, 'r') as file:
        gamma_tab = json.load(file)

    with open(R_tab_json_file_path, 'r') as file:
        R_tab = json.load(file)

    with open(T_grid_json_file_path, 'r') as file:
        T_grid = json.load(file)

    with open(phi_tab_json_file_path, 'r') as file:
        phi_tab = json.load(file)

    rho_l_tab = np.array(rho_l_tab)
    rho_v_tab = np.array(rho_v_tab)
    u_l_tab = np.array(u_l_tab)
    u_v_tab = np.array(u_v_tab)
    rho_l_tab = np.asarray(rho_l_tab)
    rho_v_tab = np.asarray(rho_v_tab)
    u_l_tab = np.asarray(u_l_tab)
    u_v_tab = np.asarray(u_v_tab)
    phi_tab = np.asarray(phi_tab)
    T_grid = np.array(T_grid)
    Tad_tab = np.array(Tad_tab)
    gamma_tab = np.array(gamma_tab)
    R_tab = np.array(R_tab)
    phi_tab = np.asarray(phi_tab)
    Tad_tab = np.asarray(Tad_tab)
    gamma_tab = np.asarray(gamma_tab)
    R_tab = np.asarray(R_tab)

    Files = {}
    Files['rho_l_tab'] = rho_l_tab
    Files['rho_v_tab'] = rho_v_tab
    Files['u_l_tab'] = u_l_tab
    Files['u_v_tab'] = u_v_tab
    Files['T_grid'] = T_grid
    Files['phi_tab'] = phi_tab
    Files['Tad_tab'] = Tad_tab
    Files['gamma_tab'] = gamma_tab
    Files['R_tab'] = R_tab

    """
    ^^^ LOAD RELEVANT DATA (DON'T DELETE)^^^
    """



    """
    vvv RUN BLOWDOWN vvv
    """

    data = blowdown_solve_set_tburn(Input, Files)

    """
    ^^^ RUN BLOWDOWN ^^^
    """



    """
    vvv ASSIGN DATA vvv
    """
    t = data[0]
    Pc = data[1]
    OF = data[2]
    phi = data[3]
    P_pressurant = data[4]
    P_injector = data[5]
    mdot_fuel_comb = data[6]
    mdot_fuel_film = data[7]
    m_fuel = data[8]
    V_fuel = data[9]
    V_pressurant = data[10]
    T_ox = data[11]
    P_ox = data[12]
    mdot_ox = data[13]
    m_ox = data[14]
    m_ox_vapor = data[15]
    m_ox_liquid = data[16]
    U_tot = data[17]
    Thrust = data[18]

    T_ox_0 = CP.PropsSI('T', 'P', Input['P_ox_0'], 'Q', 1, "N2O")
    Substance = Input['Substance']
    Ullage_ox = N2O_ullage_TVm(T_ox[0], Input['V_tank_ox'], m_ox[0])
    Ullage_fuel = V_pressurant[0]/(V_pressurant[0]+V_fuel[0])
    idx = int((1/0.001)*(-T_ox_0+217))
    V_tank_ox = Input['V_tank_ox']
    rho_l_00 = CP.PropsSI('D', 'T', T_ox_0, 'Q', 0, Substance)  # Density of liquid
    rho_v_00 = CP.PropsSI('D', 'T', T_ox_0, 'Q', 1, Substance)  # Density of vapor
    V_liquid_ox = m_ox_liquid[0] / rho_l_00
    V_vapor_ox = m_ox_vapor[0] / rho_v_00
    V_tank_ox = V_liquid_ox + V_vapor_ox
    m_ox_0 = N2O_tank_mass_TVml(T_ox_0, V_tank_ox, m_ox_liquid[0])
    rho_fuel = Input['rho_fuel']
    P_pressurant_0 = Input['P_pressurant_0']
    P_stiffness_fuel = Input['P_stiffness_fuel']
    P_ox_0 = P_ox[0]
    A_ox = Input['A_ox']
    A_inj = Input['A_inj']
    A_film = Input['A_film']
    V_tank_fuel = Input['V_tank_fuel']
    t_burn = Input['t_burn']
    """
    ^^^ ASSIGN DATA ^^^
    """



    """
    vvv PRINT DESIRED PARAMTERS vvv
    """
    init = 5
    print(f"GENERAL:")
    print(f"OF_0 = {mdot_ox[init] / mdot_fuel_comb[init]:.2f}")
    print(f"Pc_0 = {Pc[init]/6894.76:.2f} [psi]")
    print(f"mdot_fuel_0 = {mdot_fuel_comb[init]:.2f} [kg/s]")
    print(f"mdot_ox_0 = {mdot_ox[init]:.2f} [kg/s]")
    print(f"")
    print(f"FUEL:")
    print(f"Film Cooling Percentage = {mdot_fuel_comb[0]/(mdot_fuel_comb[0]+mdot_fuel_film[0]+mdot_ox[0])}")
    print(f"Film Fuel Percentage = {mdot_fuel_film[0]/(mdot_fuel_comb[0]+mdot_fuel_film[0])}%")
    print(f"Fuel Ullage = {Ullage_fuel:.2f}")
    print(f"m_fuel_0 = {m_fuel[0]:.2f} [kg]")
    print(f"V_fuel_0 = {1000*(m_fuel[0]/rho_fuel):.2f} [L]")
    print(f"P_pressurant_0 = {P_pressurant_0/6894.76:.2f} [psi]")
    print(f"P_inj_f = {P_injector[-1]/6894.76:.2f} [psi]; P_stiff_f = {(Pc[-1]+(P_stiffness_fuel*Pc[-1]))/6894.76:.2f} [psi]")
    print(f"m_fuel_f = {m_fuel[-1]:.2f} [kg]")
    print(f"")
    print(f"OX:")
    print(f"T_ox_0 = {T_ox_0:.2f} [K]")
    print(f"P_ox_0 = {P_ox_0/6894.76:.2f} [psi]")
    print(f"Ox Ullage = {Ullage_ox:.2f}")
    print(f"m_ox_liquid_0 = {m_ox_liquid[0]:.2f} [kg]")
    print(f"m_ox_liquid_0 = {m_ox_liquid[0]*2.20462:.2f} [lbm]")
    print(f"m_ox_0 = {m_ox_0:.2f} [kg]")
    print(f"m_ox_0 = {m_ox_0*2.20462:.2f} [lbm]")
    print(f"m_ox_f = {m_ox[-1]*2.20462:.2f} [lbm]")
    print(f"Tank_OF = {m_ox_liquid[0]/m_fuel[0]}")
    print(f"V_tank_ox = {V_tank_ox*1000:.2f} [L]")
    print(f"V_tank_fuel = {V_tank_fuel*1000:.2f} [L]")
    print(f"A_ox = {A_ox:.12f} [m2]")
    print(f"A_inj = {A_inj:.12f} [m2]")
    print(f"A_film = {A_film:.12f} [m2]")
    """
    ^^^ PRINT DESIRED PARAMETERS ^^^
    """



    """
    vvv SAVE DATA vvv
    """

    # json_file_ending = f"_tb{t_burn}.json"
    # save_data_json_single(data, json_file_ending)
    json_file_ending = f"_TD_HF_1.json"
    save_data_json_single(data[0], f"t{json_file_ending}")
    save_data_json_single(data[1], f"Pc{json_file_ending}")
    save_data_json_single(data[18], f"Thrust{json_file_ending}")
    save_data_json_single(data[12], f"Pox{json_file_ending}")
    save_data_json_single(data[4], f"Pfuel{json_file_ending}")
    save_data_json_single(data[14], f"mox{json_file_ending}")
    """
    ^^^ SAVE DATA ^^^
    """



    """
    vvv IT'S PLOTTING TIME vvv
    """
    generatePlot(OF, t, ["Time [s]", "OF", 1])
    generatePlot(Pc, t, ["Time [s]", "Pc [psi]", 1/6894.76])
    generatePlot(Thrust, t, ["Time [s]", "Thrust [lb]", 0.224809])

    plt.figure()
    plt.plot(t, Pc/6894.76, label='Pc')
    plt.plot(t, P_injector/6894.76, label='P_inj')
    plt.plot(t, P_pressurant/6894.76, label='P_pressurant')
    plt.plot(t, P_ox/6894.76, label='P_ox')
    plt.xlabel('time [s]')
    plt.ylabel('Pressure [psi]')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(t, m_ox, label='m_ox')
    plt.plot(t, m_ox_vapor, label='m_ox_vapor')
    plt.plot(t, m_ox_liquid, label='m_ox_liquid')
    plt.xlabel('time [s]')
    plt.ylabel('Mass [kg]')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    plt.plot(t, mdot_ox, label='mdot_ox')
    plt.plot(t, mdot_fuel_comb, label='mdot_fuel_comb')
    plt.plot(t, mdot_fuel_film, label='mdot_fuel_film')
    plt.plot(t, mdot_fuel_comb + mdot_fuel_film, label='mdot_fuel_tot')
    plt.xlabel('time [s]')
    plt.ylabel('MDot [kg/s]')
    plt.legend()
    plt.grid()
    plt.show()
    """
    ^^^ IT'S PLOTTING TIME ^^^
    """

if __name__ == "__main__":
    main()
