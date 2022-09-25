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


def calc_deltaP(mdot, Cv, T, Pin):
    """Calculates pressure drop across valve.

    Keyword arguments:
    mdot    --  
    Cv      --  
    T       --  
    Pin     --  
    """

    # calculate density in kg/m3 then convert to lbm/ft3
    dwater = CoolProp.PropsSI('D', 'T|liquid', 294, 'P', Pin, 'Water') / 16.018

    dNOx = CoolProp.PropsSI('D', 'T|liquid', T, 'P',
                            Pin, 'NitrousOxide') / 16.018

    # calculate specific gravity
    SG = dNOx / dwater

    # calculate flowrate from mdot
    Q = mdot / dNOx * 448.8

    deltaP = (Q * math.sqrt(SG) / Cv) ** 2

    return (deltaP)
