import json as json
import numpy as np


def save_data_json_single(data, json_file_ending):
    # json_file_ending = f"_tb{t_burn}.json"
    with open(f"{json_file_ending}", 'w') as f:
        json.dump(data.tolist(), f)

    return None


def save_data_json(data, json_file_ending):
    # json_file_ending = f"_tb{t_burn}.json"
    with open(f"m_ox_0_tab{json_file_ending}", 'w') as f:
        json.dump(m_ox_0_tab.tolist(), f)
    with open(f"P_ox_tab{json_file_ending}", 'w') as f:
        json.dump(P_ox_tab.tolist(), f)
    with open(f"m_ox_l_0_tab{json_file_ending}", 'w') as f:
        json.dump(m_ox_l_0_tab.tolist(), f)
    with open(f"T_ox_0_tab{json_file_ending}", 'w') as f:
        json.dump(T_ox_0.tolist(), f)
    with open(f"m_ox_min_tab{json_file_ending}", 'w') as f:
        json.dump(m_ox_min_tab.tolist(), f)

    return None


def save_data_csv(data):
    data_dict = {
        "t": t,
        "Pc": Pc,
        "OF": OF,
        "phi": phi,
        "P_pressurant": P_pressurant,
        "P_injector": P_injector,
        "mdot_fuel_comb": mdot_fuel_comb,
        "mdot_fuel_film": mdot_fuel_film,
        "m_fuel": m_fuel,
        "V_fuel": V_fuel,
        "V_pressurant": V_pressurant,
        "T_ox": T_ox,
        "P_ox": P_ox,
        "mdot_ox": mdot_ox,
        "m_ox": m_ox,
        "m_ox_vapor": m_ox_vapor,
        "m_ox_liquid": m_ox_liquid,
        "U_tot": U_tot
    }

    # Stack columns
    data_matrix = np.column_stack(list(data_dict.values()))

    # Save
    np.savetxt(
        f"BlowdownData_FilmCoolingPercentage{Finley_film_cooling_percentage:.2f}_moxliquid{m_ox_liquid[0]:.2f}_FilmUllage{Ullage_fuel:.2f}.csv",
        data_matrix,
        delimiter=",",
        header=",".join(data_dict.keys()),
        comments=""
    )

    # STEP 500
    step = 500
    data_dict = {
        "t": t[::step],
        "Pc": Pc[::step],
        "OF": OF[::step],
        "phi": phi[::step],
        "P_pressurant": P_pressurant[::step],
        "P_injector": P_injector[::step],
        "mdot_fuel_comb": mdot_fuel_comb[::step],
        "mdot_fuel_film": mdot_fuel_film[::step],
        "m_fuel": m_fuel[::step],
        "V_fuel": V_fuel[::step],
        "V_pressurant": V_pressurant[::step],
        "T_ox": T_ox[::step],
        "P_ox": P_ox[::step],
        "mdot_ox": mdot_ox[::step],
        "m_ox": m_ox[::step],
        "m_ox_vapor": m_ox_vapor[::step],
        "m_ox_liquid": m_ox_liquid[::step],
        "U_tot": U_tot[::step]
    }

    # Stack columns
    data_matrix = np.column_stack(list(data_dict.values()))

    # Save
    np.savetxt(
        f"BlowdownData_Step{step:.0f}_FilmCoolingPercentage{Finley_film_cooling_percentage:.2f}_moxliquid{m_ox_liquid[0]:.2f}_FilmUllage{Ullage_fuel:.2f}.csv",
        data_matrix,
        delimiter=",",
        header=",".join(data_dict.keys()),
        comments=""
    )

    return None
