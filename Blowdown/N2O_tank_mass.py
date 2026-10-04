import CoolProp.CoolProp as CP

def N2O_tank_mass_TVml(T, V, ml):
    Substance = "N2O"
    rho_v = CP.PropsSI('D', 'T', T, 'Q', 1, Substance)  # Density of vapor
    rho_l = CP.PropsSI('D', 'T', T, 'Q', 0, Substance)  # Density of liquid
    V_l = ml / rho_l
    V_v = V - V_l
    mv = rho_v * V_v
    m=ml+mv
    # print(f"INITIAL CONDITIONS: T={T:.2f}[K]; V={V*1000:.2f}[L]; ml={ml:.2f}[kg]")
    # print(f"MASS: m={m:.2f}[kg]; m={m*2.20462:.2f}[lbf]")

    return m

def N2O_tank_mass_PVml(P, V, ml):
    P = P
    Substance = "N2O"
    rho_v = CP.PropsSI('D', 'P', P, 'Q', 1, Substance)  # Density of vapor
    rho_l = CP.PropsSI('D', 'P', P, 'Q', 0, Substance)  # Density of liquid
    V_l = ml / rho_l
    V_v = V - V_l
    mv = rho_v * V_v
    m=ml+mv
    # print(f"INITIAL CONDITIONS: P={P/6894.76:.2f}[psi]; V={V*1000:.2f}[L]; ml={ml:.2f}[kg]")
    # print(f"MASS: m={m:.2f}[kg]; m={m*2.20462:.2f}[lbf]")
    return m

def N2O_liquid_mass_PVm(P, V, m):
    P = P
    Substance = "N2O"
    rho_v = CP.PropsSI('D', 'P', P, 'Q', 1, Substance)  # Density of vapor
    rho_l = CP.PropsSI('D', 'P', P, 'Q', 0, Substance)  # Density of liquid
    ml = m
    mv = m-ml
    Vl = ml/rho_l
    Vv = mv/rho_v
    V_temp = Vv+Vl
    err = abs((V-V_temp)/V)
    while err > 0.001:
        ml=ml-0.001
        mv=m-ml
        Vl = ml/rho_l
        Vv = mv/rho_v
        V_temp = Vv+Vl
        err = abs((V-V_temp)/V)
    U = Vv/V
    # print(f"U={Vv/V}")
    # print(f"ml={ml}")
    # print(f"mv={mv}")
    return ml
