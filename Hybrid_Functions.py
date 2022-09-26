import math
import CoolProp


def calculate_burning_area(R, r, L, rb, t):
    """Calculates burning area of fuel.

    Keyword arguments:
    R   --  OD
    r   --  ID
    L   --  length
    rb  --  burnrate
    t   --  time
    """

    if r < R:
        A = math.pi * 2 * (r + rb * t) * L
    elif r >= R:
        A = 0

    return A


def calculate_injector_oxidizer_flowrate(Cd, Ainj, dNOx, deltaP):
    """Calculates oxidizer mass flowrate through injector.

    Keyword arguments:
    Cd      --  
    Ainj    --  
    dNOx    --  
    deltaP  --  
    """

    mdot = Cd * A * math.sqrt(2 * dNOx * deltaP)

    return mdoto


def calc_mdotf(dfuel, Ab, rb):
    """Calculates fuel mass flowrate.

    Keyword arguments:
    dfuel   --  
    Ab      --  
    rb      --  
    """

    mdotf = dfuel * Ab * rb

    return mdotf


def calc_thrust(mdoto, mdotf, Isp):
    """Calculates thrust.

    Keyword arguments:
    mdoto   --  
    mdotf   --  
    Isp     --
    """

    g0 = 32.17  # ft/s^2 which is 9.81 m/s^2

    F = (mdoto + mdotf) * Isp * g0

    return F


def calculate_pressure_drop(m_dot, c_v, t, p_in, fluid):
    """Calculates the pressure drop across a valve.

    Keyword arguments:
    m_dot   --  mass flowrate (kg/s)
    c_v     --  flow coefficient
    t       --  temperature (K)
    p_in    --  inlet pressure (Pa)
    fluid   --  fluid name
    """

    WATER_DENSITY = 1000  # (kg/m^3)

    fluid_density = CoolProp.CoolProp.PropsSI('D', 'T|liquid', t, 'P',
                                     p_in, fluid)  # Fluid density (kg/m^3)

    fluid_specific_gravity = fluid_density / \
        WATER_DENSITY  # Fluid specific gravity

    q = m_dot / fluid_density  # Volume flowrate (m^3/s)

    q_gal_per_min = q * 15850  # Volume flowrate (gal/min)

    delta_p_psi = (q_gal_per_min * math.sqrt(fluid_specific_gravity) /
                   c_v) ** 2  # Pressure drop (psi)

    delta_p = delta_p_psi * 6895  # Pressure drop (Pa)

    return delta_p
