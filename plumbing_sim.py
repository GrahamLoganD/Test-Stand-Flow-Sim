import math
import CoolProp
import hybrid_functions

fluid = 'NitrousOxide'

atmospheric_pressure = 101325  # Surrounding atmospheric pressure (Pa)
bottle_weight_lbf = 10  # Weight of fluid in the bottle (lbf)
bottle_pressure_psi = 750  # Pressure inside the bottle (psig)
atmospheric_temperature_F = 68  # Surrounding atmospheric temperature (F)

# Surrounding atmospheric temperature (K)
atmospheric_temperature = (atmospheric_temperature_F - 32) * 5 / 9 + 273.15

# Mass of fluid in the bottle (kg)
bottle_mass = bottle_weight_lbf * 0.45359237

bottle_pressure = bottle_pressure_psi * 6895  # Pressure inside the bottle (Pa)

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
reqdeltaP = (mdot / (Cd * Ainj)) ** 2 / (2 * (CoolProp.CoolProp.PropsSI('D',
                                                                        'T', T, 'P', P3 * 6895, fluid) / 16.018 * 32.2))
Pc = P4 - reqdeltaP
# print results
print(P1, P2, P3, P4, Pc)
print(reqdeltaP)
