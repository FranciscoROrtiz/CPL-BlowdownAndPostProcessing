# filename: Cd_and_dP_avg_calc.py

import math as math
import numpy as np

def Cd_and_dP_avg_calc(t, p1, p2, m_discharged, rho, A):
    """
    Calculate the Cd

    Assume SPI model across injector orifices in order to calculate Cd from
    upstream pressure, downstream pressure, and mass discharged test data

    Paramters
    ---------
    t : numpy.ndarray
        Time array -> [s]
    p1 : numpy.ndarray
        Injector upstream pressure array corresponding to ``t`` -> [Pa]
    p2 : numpy.ndarray
        Injector downstream pressure array corresponding to ``t`` -> [Pa]
    m_discharged : float
        Total mass discharged over time interval -> [kg]
    rho : numpy.ndarray
        Density of fluid corresponding to ``p1`` -> [kg/m3]
    A : float
        Area of injector holes -> [m2]

    Returns
    -------
    Cd : float
        Discharge coefficient
    dP_avg : float
        Average p1-p2 with respect to time

    """
    
    dP = p1 - p2
    integral = np.trapz(np.sqrt(rho*dP), t)
    m_discharged_ideal = math.sqrt(2)*A*integral

    Cd = m_discharged / m_discharged_ideal
    dP_avg = (1/(t[-1]-t[0]))*np.trapz(dP, t)

    return Cd, dP_avg
