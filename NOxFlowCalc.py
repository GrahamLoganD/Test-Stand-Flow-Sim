import math
import CoolProp

# Initialize Independent Variables
T = 294  # 70f=294K
P = 800*6895  # psi to pa
Cv = 0.4  # from valve specs
DP = 1*6895  # psi to pa

# Calculate Dependent Variables
# calculate density in kg/m3 then convert to lbm/ft3
dwater = CoolProp.PropsSI('D', 'T|liquid', T, 'P', P, 'Water')/16.018
dNOx = CoolProp.PropsSI('D', 'T', T, 'P', P, 'NitrousOxide')/16.018
# calculate specific gravity
SG = dNOx/dwater
# calcualte volume flowrate in GPM
Qgpm = Cv/math.sqrt(SG/DP)
# convert GPM to ft3/s
Q = Qgpm*0.002228
mdot = dNOx*Q

# Display Important Variables
print("SG:", SG)
print("mdot:", mdot, "lbs/s")
