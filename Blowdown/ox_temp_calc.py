import numpy as np

def ox_temp_calc(U_tot, m_ox, V_tank, rho_l_tab, rho_v_tab, u_l_tab, u_v_tab, T_grid):
    x = (U_tot / m_ox - u_l_tab) / (u_v_tab - u_l_tab)
    V_tank_temp = m_ox * (((1 - x) / rho_l_tab) + (x / rho_v_tab))

    V_tank_temp = np.asarray(V_tank_temp)
    # print(x)
    # # print(V_tank_temp)
    # plt.figure()
    # plt.plot(T_grid, V_tank_temp*1000)
    # plt.xlabel('T_grid [K]')
    # plt.ylabel('V_tank_temp [L]')
    # plt.grid()
    # plt.show()

    idx = (np.abs(V_tank_temp - V_tank)).argmin()
    T_new = T_grid[idx]
    x_new = x[idx]
    # print(f"T_new = {T_new}")
    return T_new, x_new
