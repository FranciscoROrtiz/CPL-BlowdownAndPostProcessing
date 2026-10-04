import CoolProp.CoolProp as CP

def N2O_ullage_TVm(T, V, m):
    Substance = "N2O"
    rho_v = CP.PropsSI('D', 'T', T, 'Q', 1, Substance)  # Density of vapor
    rho_l = CP.PropsSI('D', 'T', T, 'Q', 0, Substance)  # Density of liquid
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
    U=Vv/V
    # print(f"U={U}")
    # print(f"ml={ml}")
    # print(f"mv={mv}")
    return U

def N2O_ullage_PVm(P, V, m):
    P = P*6894.76
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
    return U
