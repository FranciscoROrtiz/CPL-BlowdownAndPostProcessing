import math as math
from def_combustion_gas_new import def_combustion_gas_new

def calc_chamber_pressure(phi, mdot_total, At, phi_tab, Tad_tab, gamma_tab, R_tab, n_cstar):
    '''
    Inputs: phi = Equivalence Ratio, pc = Chamber Pressure [Pa], mdot_total = Total Mass Flow Rate (fuel + ox)
    Outputs: pc = chamber pressure
    HOW TO SOLVE FOR CHAMBER PRESSURE
    c*_t = sqrt(R*gc*Tc) / GAMMA
    c*_e = At * pc * gc / mdot

    1. c*_t = c*_e --> sqrt(R*gc*Tc) / GAMMA = At * pc * gc / mdot
    2. Solve for pc --> pc = (sqrt(R*gc*Tc) / GAMMA) * (mdot / At * gc)
    '''
    phi0 = 6.59 / 3      # Design OF of 3 at t = 0s
    p0 = 2.0684 * 10**6     # Design pc of 300 psi at t = 0s
    # Define Gas at new Equivalence Ratio, Chamber Pressure
    Tc, gamma, R = def_combustion_gas_new(phi, phi_tab, Tad_tab, gamma_tab, R_tab)     # Chamber Temperature, Ratio of Specific Heats, Mixture Gas Constant
    GAMMA = math.sqrt(gamma * ((gamma + 1) / 2)**( -1 * (gamma + 1) / (gamma - 1) ) )       # Flow Function
    gc = 1       # %pretty sure g represents an imperiall unit conversion --> default to 1 here bc we use SI units

    cstar_t = math.sqrt(R*gc*Tc) / GAMMA        # Theoretical c* formula; c*_t(phi_i-1, pc_i-1)
    pc_new = cstar_t * ( mdot_total / (At * gc) ) * (n_cstar)   # New Chamber Pressure

    return pc_new
