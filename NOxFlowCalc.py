import CoolProp
import math
from math import sqrt
from CoolProp.CoolProp import PropsSI

# Initialize Independent Variables
T = 294  # 70f=294K
P = 800*6895  # psi to pa
Cv = 0.4  # from valve specs
DP = 1*6895  # psi to pa

# Calculate Dependent Variables
# calculate density in kg/m3 then convert to lbm/ft3
dwater = PropsSI('D', 'T|liquid', T, 'P', P, 'Water')/16.018
dNOx = PropsSI('D', 'T', T, 'P', P, 'NitrousOxide')/16.018
# calculate specific gravity
SG = dNOx/dwater
# calcualte volume flowrate in GPM
Qgpm = Cv/sqrt(SG/DP)
# convert GPM to ft3/s
Q = Qgpm*0.002228
mdot = dNOx*Q

# Display Important Variables
print("SG:", SG)
print("mdot:", mdot, "lbs/s")
