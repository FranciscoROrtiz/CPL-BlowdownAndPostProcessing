def mdot_SPI(Cd, A, rho, dP):
    # print(f"Cd={Cd}")
    # print(f"A={A:.15f}")
    # print(f"rho={rho}")
    # print(f"dP={dP}")
    mdot = Cd*A*(2*rho*dP)**0.5
    return mdot
