import math
import CoolProp
import hybrid_functions

FLUID = 'NitrousOxide'

ATMOSPHERIC_PRESSURE = 101325  # Surrounding atmospheric pressure (Pa)
ATMOSPHERIC_TEMPERATURE_F = 68  # Surrounding atmospheric temperature (F)
INITIAL_BOTTLE_WEIGHT_LBF = 10  # Weight of fluid in the bottle (lbf)
INITIAL_BOTTLE_PRESSURE_PSI = 750  # Pressure inside the bottle (psig)

# Flow coefficients
flow_coefficient_tank_valve = 0.41 # https://catalog.circlevalve.com/item/check-valves/200-series-check-valves-0-to-3000-psig/cs-259a-1pp
flow_coefficient_check_valve = 1.6

loss_coefficient_tee = 0.9  # Line flow, threaded
loss_coefficient_ball_valve = 0.05  # Fully open

pipe_inner_diameter = 0.00707  # Inner diamater of the pipe (m)

# Surrounding atmospheric temperature (K)
atmospheric_temperature = (ATMOSPHERIC_TEMPERATURE_F - 32) * 5 / 9 + 273.15

# Mass of fluid in the bottle (kg)
bottle_mass = INITIAL_BOTTLE_WEIGHT_LBF * 0.45359237

# Temperature inside the bottle (K)
bottle_temperature = atmospheric_temperature
bottle_pressure = INITIAL_BOTTLE_PRESSURE_PSI * 6895  # Pressure inside the bottle (Pa)

print("The bottle phase is", CoolProp.CoolProp.PhaseSI('P', bottle_pressure, 'T', bottle_temperature, FLUID),
      "and the exit phase is", CoolProp.CoolProp.PhaseSI('P', ATMOSPHERIC_PRESSURE, 'T', atmospheric_temperature, FLUID))

# Isentropic incompressible flow assumption
isentropic_incompressible_flow_velocity = math.sqrt(2 / CoolProp.CoolProp.PropsSI('D', 'T', atmospheric_temperature, 'P',
                                                                                  (bottle_pressure - ATMOSPHERIC_PRESSURE) / 2, FLUID) * (bottle_pressure - ATMOSPHERIC_PRESSURE))
isentropic_incompressible_mass_flowrate = math.pi * (pipe_inner_diameter / 2) ** 2 * isentropic_incompressible_flow_velocity * CoolProp.CoolProp.PropsSI('D', 'T', atmospheric_temperature, 'P',
                                                                                                                                                         (bottle_pressure - ATMOSPHERIC_PRESSURE) / 2, FLUID)
print("The isentropic incompressible flowrate is",
      round(isentropic_incompressible_mass_flowrate, 3) / 0.45359237, "lbm/s")

# Incompressible flow assumption
def calculate_total_pressure_drop(mass_flowrate_guess):
    pressure_drop_tank_valve = hybrid_functions.calculate_pressure_drop(
        mass_flowrate_guess, flow_coefficient_tank_valve, bottle_temperature, bottle_pressure, FLUID)

    tee_1_pressure = bottle_pressure - pressure_drop_tank_valve
    tee_1_pressure_drop = 0

    check_valve_pressure = tee_1_pressure - tee_1_pressure_drop
    check_valve_pressure_drop = hybrid_functions.calculate_pressure_drop(
        mass_flowrate_guess, flow_coefficient_check_valve, bottle_temperature, bottle_pressure, FLUID)

    tee_2_pressure = check_valve_pressure - check_valve_pressure_drop
    tee_2_pressure_drop = 0

    ball_valve_pressure = tee_2_pressure - tee_2_pressure_drop
    ball_valve_pressure_drop = 0

    tee_3_pressure = ball_valve_pressure - ball_valve_pressure_drop
    tee_3_pressure_drop = 0

    exit_pressure = tee_3_pressure - tee_3_pressure_drop

    return bottle_pressure - exit_pressure


mass_flowrate_upper = isentropic_incompressible_mass_flowrate
mass_flowrate_lower = 0

pressure_drop_upper = calculate_total_pressure_drop(mass_flowrate_upper)
pressure_drop_lower = calculate_total_pressure_drop(mass_flowrate_lower)

real_pressure_drop = bottle_pressure - ATMOSPHERIC_PRESSURE

error_upper = pressure_drop_upper - real_pressure_drop
error_lower = pressure_drop_lower - real_pressure_drop

error_guess = 1

# Bisection root-finding method
while abs(error_guess) > (10 ** -3):
    mass_flowrate_guess = (mass_flowrate_lower + mass_flowrate_upper) / 2
    pressure_drop_guess = calculate_total_pressure_drop(mass_flowrate_guess)
    error_guess = pressure_drop_guess - real_pressure_drop

    if error_lower * error_guess < 0:
        mass_flowrate_upper = mass_flowrate_guess
    elif error_lower * error_guess > 0:
        mass_flowrate_lower = mass_flowrate_guess

incompressible_mass_flowrate = mass_flowrate_guess

print("The incompressible flowrate is",
      round(incompressible_mass_flowrate, 3) / 0.45359237, "lbm/s")
print("Pressure compare", round(pressure_drop_guess, 3) / (10 ** 3), "kPa to",
      round(bottle_pressure - ATMOSPHERIC_PRESSURE, 3) / (10 ** 3), "kPa")

""""
# Calculate Pressure Drop Out of Tank S/O
deltaP_tank = hybrid_functions.calculate_pressure_drop(
    mdot, .41, T, Ptank * 6895, FLUID)
P1 = Ptank - deltaP_tank
# Calculate Pressure Drop Across HOV-01
deltaP_HOV01 = hybrid_functions.calculate_pressure_drop(
    mdot, 6, T, P1 * 6895, FLUID)
P2 = P1 - deltaP_HOV01
# Calculate Pressure Drop Across CV-01
deltaP_CV01 = hybrid_functions.calculate_pressure_drop(
    mdot, 3.1, T, P2 * 6895, FLUID)
P3 = P2 - deltaP_CV01
# Calculate Pressure Drop Across HOV-02
deltaP_HOV02 = hybrid_functions.calculate_pressure_drop(
    mdot, 6, T, P3 * 6895, FLUID)
P4 = P3 - deltaP_HOV02

Cd = 0.7
rports = 0.03
numports = 25
Ainj = math.pi * rports ** 2 * numports
reqdeltaP = (mdot / (Cd * Ainj)) ** 2 / (2 * (CoolProp.CoolProp.PropsSI('D',
                                                                        'T', T, 'P', P3 * 6895, FLUID) / 16.018 * 32.2))
Pc = P4 - reqdeltaP
# print results
print(P1, P2, P3, P4, Pc)
print(reqdeltaP)
"""
