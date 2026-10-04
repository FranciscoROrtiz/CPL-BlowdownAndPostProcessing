def calculate_injector_Area(G, mdot, Cd):
    Area = mdot / (G * Cd)

    return Area
