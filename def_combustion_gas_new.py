import numpy as np

def def_combustion_gas_new(phi, phi_tab, Tad_tab, gamma_tab, R_tab):
    '''
    Inputs: phi = Equivalence Ratio, pc = Chamber Pressure [Pa]
    Outputs: Tad = Adiabatic Flame Temperature [K], gamma_mix = Ratio of specific heats of mixture, R_mix = Mixture Gas Constant
    :return:
    '''

    Tad = np.interp(phi, phi_tab, Tad_tab)
    gamma_mix = np.interp(phi, phi_tab, gamma_tab)
    R_mix = np.interp(phi, phi_tab, R_tab)
    return Tad, gamma_mix, R_mix
