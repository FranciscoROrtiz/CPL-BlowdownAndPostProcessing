from scipy.optimize import brentq
from def_combustion_gas_new import def_combustion_gas_new
from area_mach_residual import area_mach_residual

def thrust(Pc, phi, mdot_total, Ae, A_th, phi_tab, Tad_tab, gamma_tab, R_tab):
    Tad, gamma, R = def_combustion_gas_new(phi, phi_tab, Tad_tab, gamma_tab, R_tab)     # Chamber Temperature, Ratio of Specific Heats, Mixture Gas Constant
    # print(Tad)
    Pa = 101325
    # Pa = 61000
    target_AR = Ae/A_th
    # target_AR = 3.5977363508
    Ma = brentq(area_mach_residual, 1.0, 20.0, args=(target_AR, gamma))
    T_ratio = (1+(((gamma-1)/2)*Ma))
    P_ratio = (1+(((gamma-1)/2)*Ma))**(gamma/(gamma-1))
    Te = Tad/T_ratio
    Pe = Pc/P_ratio
    rhoe = Pe/(R*Te)

    Ve = mdot_total/(rhoe*Ae)

    Thrust = mdot_total*Ve + (Pe-Pa)*Ae

    return Thrust
