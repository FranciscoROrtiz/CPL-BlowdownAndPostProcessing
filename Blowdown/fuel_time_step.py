from mdot_SPI import mdot_SPI
from AdiabaticBlowdown import AdiabaticBlowdown

def fuel_time_step(P_loss_fuel, Cd_fuel, A_inj, A_film, rho_fuel, dt, P_pressurant_curr, m_fuel_curr, V_fuel_curr, V_pressurant_curr, gamma, Pc_curr):
    dP = P_pressurant_curr - P_loss_fuel - Pc_curr
    # if dP<=0:
    #     dP = 101325
    # print(f"Pc = {Pc_curr/6894.76:.2f}")
    # print(f"P_pressurant_curr = {P_pressurant_curr/6894.76}")
    # print(f"P_pressurant_curr={P_pressurant_curr/6894.76:.2f}")
    # print(f"Pc_curr = {Pc_curr/6894.76:.2f}")
    # print(1)
    mdot_fuel_comb_new = mdot_SPI(Cd_fuel, A_inj, rho_fuel, dP)
    # print(2)
    mdot_fuel_film_new = mdot_SPI(Cd_fuel, A_film, rho_fuel, dP)

    dm = (mdot_fuel_comb_new + mdot_fuel_film_new) * dt
    m_fuel_new = m_fuel_curr - dm

    V_fuel_new = m_fuel_new / rho_fuel

    dV = V_fuel_curr - V_fuel_new
    V_pressurant_new = V_pressurant_curr + dV

    P_pressurant_new = AdiabaticBlowdown(V_pressurant_curr, V_pressurant_new, P_pressurant_curr, gamma)
    P_pressurant_new = P_pressurant_new - (5*6894.76*dt)
    P_injector_new = P_pressurant_new - P_loss_fuel

    return P_pressurant_new, P_injector_new, mdot_fuel_comb_new, mdot_fuel_film_new, m_fuel_new, V_fuel_new, V_pressurant_new
