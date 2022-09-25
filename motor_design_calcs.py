import math
import CoolProp
import Hybrid_Functions

Pc = 400  # chamber pressure in psi
Cf = 1  # thrust coefficient
F = 100  # lbf
mr = 7.67
Isp = 186.29
g0 = 32.17
cstar = 0.93 * Isp * g0 / Cf

At = F / (Cf * Pc)
mdot = g0 * Pc * At / cstar
mdotf = mdot / (mr + 1)
mdoto = mdot - mdotf

print(mdoto, "lbs/s N2O")
print(mdotf, "lbs/s HTPB needed")

R = 2.5  # outer grain diamter
r = .1  # core diamter
L = 6  # grain length
T = 294
G0 = mdoto / (math.pi * r ** 2) * 703.1  # lbm/(s*in2) to kg/(s*m2)
rb = ((9.3368 * 10 ** -8) * G0 ** 1.6386) * 3.281  # m/s to ft/s
t = 0
Ab = Hybrid_Functions.calc_Ab(R, r, L, rb, t) / 144  # in2 to ft2
provided_mdotf = Hybrid_Functions.calc_mdotf(
    103.57, Ab, rb)  # google says 57.43 lb/ft3
print(provided_mdotf, "provided mdotf")
