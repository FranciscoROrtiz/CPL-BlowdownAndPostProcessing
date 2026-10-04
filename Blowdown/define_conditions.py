import CoolProp.CoolProp as CP

def define_conditions(p_1, p_2):
    # Define the fluid and conditions
    # print(p_2)
    fluid = 'N2O'  # Nitrous Oxide
    fluid_properties = {
        'P_1': p_1,  # Pressure in Pascals
        'P_2': p_2   # Pressure in Pascals
    }

    # Instead of assuming T_1, get it from saturation at P_1
    fluid_properties['T_1'] = CP.PropsSI('T', 'P', fluid_properties['P_1'], 'Q', 0, fluid)

    # Use saturated liquid entropy at upstream pressure
    fluid_properties['entropy'] = CP.PropsSI('S', 'P', fluid_properties['P_1'], 'Q', 0, fluid)
    # print(fluid_properties['P_2'])
    fluid_properties['density_1'] = CP.PropsSI('D', 'S', fluid_properties['entropy'], 'P', fluid_properties['P_1'], fluid)
    # fluid_properties['density_2'] = CP.PropsSI('D', 'S', fluid_properties['entropy'], 'P', fluid_properties['P_2'], fluid)
    # fluid_properties['enthalpy_1'] = CP.PropsSI('H', 'S', fluid_properties['entropy'], 'P', fluid_properties['P_1'], fluid)
    # fluid_properties['enthalpy_2'] = CP.PropsSI('H', 'S', fluid_properties['entropy'], 'P', fluid_properties['P_2'], fluid)
    # fluid_properties['T_2'] = CP.PropsSI('T', 'S', fluid_properties['entropy'], 'P', fluid_properties['P_2'], fluid)

    # Vapor pressure at T_1 (for use in SPI or Dyer model if needed)
    fluid_properties['P_v'] = CP.PropsSI('P', 'T', fluid_properties['T_1'], 'Q', 0, fluid)

    return fluid_properties
