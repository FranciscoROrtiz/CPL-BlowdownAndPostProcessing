def calculate_mass_flux_rate(N2O):
    """
    :param N2O: Dictionary of N2O properties
    :return: Critical Mass Flux Rate of Nitrous in kg/(m^2*s)
    """
    # k = ((P_1-P_2)/(P_v-P_2))**0.5
    # G_SPI = (2*density*(P_1-P_2))**0.5  # Calculate the SPI critical mass flux rate
    # G_HEM = density*((2*(enthalpy_1-enthalpy_2))**0.5)  # Calculate the HEM critical mass flux rate

    # k = ((N2O['P_1'] - N2O['P_2']) / (N2O['P_v'] - N2O['P_2'])) ** 0.5
    G_SPI = (2*N2O['density_1']*(N2O['P_1']-N2O['P_2']))**0.5  # Calculate the SPI critical mass flux rate
    # G_HEM = N2O['density_2']*((2*(N2O['enthalpy_1']-N2O['enthalpy_2']))**0.5)  # Calculate the HEM critical mass flux rate
    # G_DYER = ((k/(1+k))*G_SPI + (1/(1+k))*G_HEM)  # Calculate the DYER critical mass flux rate

    return G_SPI
