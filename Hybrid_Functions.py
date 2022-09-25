import math
import CoolProp


def calc_Ab(R, r, L, rb, t):
    # calcualtes burning area of fuel
    #R=OD, r=ID, L=length, rb=burnrate, t=time
    if r < R:
        A = math.pi * 2 * (r + rb * t) * L
    elif r >= R:
        A = 0
    return A


def calc_mdoto(Cd, Ainj, dNOx, deltaP):
    # calculates oxidizer mass flowrate through injecttor
    mdot = Cd * A * math.sqrt(2 * dNOx * deltaP)
    return mdoto


def calc_mdotf(dfuel, Ab, rb):
    # calculates fuel mass flowrate
    mdotf = dfuel * Ab * rb
    return mdotf


def calc_thrust(mdoto, mdotf, Isp):
    # calculates thrust
    g0 = 32.17  # ft/s^2 which is 9.81 m/s^2
    F = (mdoto + mdotf) * Isp * g0
    return F


def calc_deltaP(mdot, Cv, T, Pin):
    # calculates pressure drop across valve
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
