import math
import CoolProp
import hybrid_functions

fluid = 'NitrousOxide'
mdot = .5
Ptank = 750
T = 294
# Calculate Tank Pressure
mbot = 10  # lbm of N2O in full bottle
Ttank = T

# Calculate Pressure Drop Out of Tank S/O
deltaP_tank = hybrid_functions.calculate_pressure_drop(
    mdot, .41, T, Ptank * 6895, fluid)
P1 = Ptank - deltaP_tank
# Calculate Pressure Drop Across HOV-01
deltaP_HOV01 = hybrid_functions.calculate_pressure_drop(
    mdot, 6, T, P1 * 6895, fluid)
P2 = P1 - deltaP_HOV01
# Calculate Pressure Drop Across CV-01
deltaP_CV01 = hybrid_functions.calculate_pressure_drop(
    mdot, 3.1, T, P2 * 6895, fluid)
P3 = P2 - deltaP_CV01
# Calculate Pressure Drop Across HOV-02
deltaP_HOV02 = hybrid_functions.calculate_pressure_drop(
    mdot, 6, T, P3 * 6895, fluid)
P4 = P3 - deltaP_HOV02

Cd = 0.7
rports = 0.03
numports = 25
Ainj = math.pi * rports ** 2 * numports
reqdeltaP = (mdot / (Cd * Ainj)) ** 2 / (2 * (CoolProp.PropsSI('D',
                                                               'T', T, 'P', P3 * 6895, fluid) / 16.018 * 32.2))
Pc = P4 - reqdeltaP
# print results
print(P1, P2, P3, P4, Pc)
print(reqdeltaP)
