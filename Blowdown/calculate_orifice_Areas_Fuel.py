def calculate_orifice_Areas_Fuel(mdot_fuel, Cd, p1, p2, rho, film_cooling_percentage):
    """
    :param p1: Injector Inlet Pressure, p2: Chamber Pressure,
    rho: Density of Fuel, Cd: Discharge coefficient, mdot_fuel: mass flow rate of fuel
    :return: A_injector: total Area for fuel injection, A_film, total area for film cooling
    """

    # film_cooling_percentage = 0.2

    # G_SPI = (2*density*(P_1-P_2))**0.5  # Calculate the SPI critical mass flux rate

    G_SPI = (2*rho*(p1-p2))**0.5  # Calculate the SPI critical mass flux rate
    A_injector = mdot_fuel * (1 - film_cooling_percentage) / (G_SPI * Cd)
    A_film = mdot_fuel * film_cooling_percentage / (G_SPI * Cd)


    return A_injector, A_film
