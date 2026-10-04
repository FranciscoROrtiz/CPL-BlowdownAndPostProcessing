import CoolProp.CoolProp as CP

from define_conditions import define_conditions
from calculate_mass_flux_rate import calculate_mass_flux_rate
from ox_temp_calc import ox_temp_calc

def ox_time_step(P_loss_ox, Cd_ox, A_ox, Substance, V_tank_ox, dt, P_ox_curr, m_ox_curr, T_ox_curr, U_tot_curr, Pc_curr, rho_l_tab, rho_v_tab, u_l_tab, u_v_tab, T_grid):
    N2O = define_conditions((P_ox_curr - P_loss_ox), Pc_curr)
    G = calculate_mass_flux_rate(N2O)
    mdot_ox_new = G*Cd_ox*A_ox
    # Choked
    gamma = CP.PropsSI('CPMASS', 'T', T_ox_curr, 'Q', 1, Substance) / CP.PropsSI('CVMASS', 'T', T_ox_curr, 'Q', 1, Substance)
    P_ratio = P_ox_curr*((gamma+1)/2)**(-gamma/(gamma-1))
    rho_crit = CP.PropsSI('D', 'P', P_ratio, 'Q', 1, Substance)
    de = 0.055*0.0254
    Ae = (3.14/4)*(de)**2
    T_crit = CP.PropsSI('T', 'P', P_ratio, 'Q', 1, Substance)
    R = 188.91
    c = (gamma*R*T_crit)**0.5
    mdot_ox_choked = rho_crit*c*Ae
    # print(f"mdot_ox_new = {mdot_ox_new}")
    # print(f"mdot_ox_choked = {mdot_ox_choked}")
    # mdot_ox_choked = 0
    # Choked
    # print(f"mdot_ox_new = {mdot_ox_new}")
    dm = -mdot_ox_new * dt + (-mdot_ox_choked * dt)
    # print(f"dm = {dm}")
    m_ox_new = m_ox_curr + dm
    # print(f"m_ox_curr = {m_ox_curr}")
    # print(f"m_ox_new = {m_ox_new}")
    h_out = CP.PropsSI('H', 'T', T_ox_curr, 'Q', 0, Substance)  # Enthalpy of liquid
    h_out_choked = CP.PropsSI('H', 'T', T_ox_curr, 'Q', 1, Substance)
    dU_tot = (-mdot_ox_new * h_out) * dt + (-mdot_ox_choked * h_out_choked * dt)
    U_tot_new = U_tot_curr + dU_tot

    T_ox_new, x_new = ox_temp_calc(U_tot_new, m_ox_new, V_tank_ox, rho_l_tab, rho_v_tab, u_l_tab, u_v_tab, T_grid)


    P_ox_new = CP.PropsSI('P', 'T', T_ox_new, 'Q', 1, Substance)
    m_ox_vapor_new = x_new * m_ox_new
    m_ox_liquid_new = (1 - x_new) * m_ox_new

    return T_ox_new, P_ox_new, mdot_ox_new, m_ox_new, m_ox_vapor_new, m_ox_liquid_new, U_tot_new
