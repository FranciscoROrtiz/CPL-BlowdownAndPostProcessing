# filename: film_cooling_percentage_conversion.py



def film_cooling_percent_fuel_to_film_cooling_percent_total(mdot_ox, mdot_fuel_comb, film_cooling_percent_fuel):
    """
    Convert mdot_fuel_film/mdot_fuel --> mdot_fuel_film / mdot_total

    Paramters
    ---------
    mdot_ox : float
        Oxidizer mass flow rate -> [kg/s]
    mdot_fuel_comb : float
        Fuel mass flow rate that undergoes combustion -> [kg/s]
    film_cooling_percent_fuel : float
        Ratio of fuel mass flow rate used for film cooling and total fuel mass
        flow rate

    Returns
    -------
    float
        Ratio of fuel mass flow rate used for film cooling and total mass
        flow rate

    """

    film_cooling_percent_total = (film_cooling_percent_fuel * mdot_fuel_comb) / (mdot_ox + mdot_fuel_comb - (film_cooling_percent_fuel*mdot_ox))

    return film_cooling_percent_total

def film_cooling_percent_total_to_film_cooling_percent_fuel(mdot_ox, mdot_fuel_comb, film_cooling_percent_total):
    """
    Convert mdot_fuel_film / mdot_total --> mdot_fuel_film / mdot_fuel

    Paramters
    ---------
    mdot_ox : float
        Oxidizer mass flow rate -> [kg/s]
    mdot_fuel_comb : float
        Fuel mass flow rate that undergoes combustion -> [kg/s]
    film_cooling_percent_total : float
        Ratio of fuel mass flow rate used for film cooling and total mass
        flow rate

    Returns
    -------
    float
        Ratio of fuel mass flow rate used for film cooling and total fuel mass
        flow rate

    """

    film_cooling_percent_total = 1 - ( ( (1-film_cooling_percent_total) * mdot_fuel_comb) / ( (film_cooling_percent_total*mdot_ox) + mdot_fuel_comb) )

    return film_cooling_percent_total
