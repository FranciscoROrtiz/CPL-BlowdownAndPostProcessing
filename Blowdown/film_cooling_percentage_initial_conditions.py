# filename: film_cooling_percentage.py

from Cd_and_dP_avg_calc import Cd_and_dP_avg_calc

def film_cooling_percentage_initial_conditions(Cd_ox, A_ox, P1_ox, Cd_fuel, A_fuel, P1_fuel):
    """
    Calculate film cooling percentage for a determined set of oxidizer, fuel combusted, and film cooling mass flow rates

    Paramters
    ---------
    mdot_ox : float
        Oxidizer mass flow rate -> [kg/s]
    mdot_fuel_comb : float
        Fuel mass flow rate that undergoes combustion -> [kg/s]
    mdot_fuel_film : float
        Fuel mass flow rate that is used for film cooling -> [kg/s]

    Returns
    -------
    float
        film cooling percentage

    """

    film_cooling_percentage = mdot_fuel_film / (mdot_ox + mdot_fuel_comb + mdot_fuel_film)

    return film_cooling_percentage
