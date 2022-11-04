import CoolProp
import math
import scipy

FLUID = 'NitrousOxide'

INITIAL_BOTTLE_WEIGHT_LBF = 10  # Weight of fluid in the bottle (lbf)
INITIAL_BOTTLE_GAUGE_PRESSURE_PSI = 750  # Pressure inside the bottle (psig)

ATMOSPHERIC_PRESSURE = 101325  # Surrounding atmospheric pressure (Pa)
ATMOSPHERIC_TEMPERATURE_F = 68  # Surrounding atmospheric temperature (F)

PIPE_INNER_DIAMETER = 0.00707  # Inner diamater of the pipe (m)

# Flow coefficients
# https://catalog.circlevalve.com/item/check-valves/200-series-check-valves-0-to-3000-psig/cs-259a-1pp
TANK_VALVE_FLOW_COEFFICIENT = 0.41
CHECK_VALVE_FLOW_COEFFICIENT = 1.6

# Loss coefficients
TEE_LOSS_COEFFICIENT = 0.9  # Line flow, threaded
BALL_VALVE_LOSS_COEFFICIENT = 0.05  # Fully open


def convert_f_to_k(temperature_f):
    """Converts temperature in degrees Fahrenheit into kelvins.

    Keyword arguments:
    temperature_f   --  temperature (f)
    """

    return (temperature_f + 459.67) / 1.8


def convert_lb_to_kg(mass_lb):
    """Converts mass in pounds into kilograms.

    Keyword arguments:
    mass_lb --  mass (lb)
    """

    return 0.45359237 * mass_lb


def convert_psi_to_pa(pressure_psi):
    """Converts pressure in pounds per square inch into pascals.

    Keyword arguments:
    pressure_psi    --  pressure (psi)
    """

    return 6.894757e3 * pressure_psi


def convert_gauge_pressure_to_absolute(pressure, atmospheric_pressure):
    """Converts gauge pressure in pascals into absolute pressure.

    Keyword arguments:
    pressure                --  pressure (Pa)
    atmospheric_pressure    --  atmospheric pressure (Pa)
    """

    return pressure + atmospheric_pressure


def calculate_pipe_area(inner_diameter):
    """Calculates the cross-sectional area of the pipe.

    Keyword arguments:
    inner_diameter  --  inner diameter of the pipe (m)
    """

    return math.pi * (inner_diameter / 2) ** 2


def bernoulli_equation_exit_velocity(pressure_in, velocity_in, pressure_out, density):
    """Calculates the exit flow velocity in m/s using the Bernoulli equation, assuming the fluid is incompressible and there are no friction losses.

    Keyword arguments:
    pressure_in --  entry pressure (Pa)
    velocity_in --  entry velocity (m/s)
    pressure_in --  exit pressure (Pa)
    density     --  density of the fluid (kg/m^3)
    """

    # The energy per unit mass of the fluid at the entry point
    energy_in = pressure_in + 1 / 2 * density * velocity_in ** 2

    # The velocity of the fluid at the exit point
    return ((energy_in - pressure_out) * 2 / density) ** 0.5


def calculate_mass_flowrate_from_velocity(velocity, density, pipe_area):
    """Calculates the mass flowrate of a fluid.

    Keyword arguments:
    velocity    --  flow velocity of the fluid (m/s)
    density     --  density of the fluid (kg/m^3)
    pipe_area   --  cross-sectional area of the fluid flow (m^2)
    """

    return pipe_area * velocity * density


def convert_kg_per_s_to_lbm_per_s(mass_flowrate):
    """Converts mass flowrate in kg/s into lbm/s.

    Keyword arguments:
    mass_flowrate   --  mass flowrate of the fluid (kg/s)
    """

    return mass_flowrate / 0.45359237


def calculate_specific_gravity(density):
    """Calculates the specific gravity of a fluid.

    Keyword arguments:
    density --  fluid density (kg/m^3)
    """

    WATER_DENSITY = 1000  # (kg/m^3)

    return density / \
        WATER_DENSITY


def convert_cubic_meters_per_second_to_gpm(volume_flowrate):
    """Converts volume flowrate in m^3/s into gal/min.

    Keyword arguments:
    volume_flowrate --  volume flowrate (m^3/s)
    """

    return 1.585e4 * volume_flowrate


def calculate_flow_coefficient_pressure_drop(flow_coefficient, volume_flowrate, specific_gravity):
    """Calculates the pressure drop across a valve using the flow coefficient.

    Keyword arguments:
    flow_coefficient    --  flow coefficient of the valve
    volume_flowrate     --  rate of flow of the fluid (m^3/s)
    specific gravity    --  specific gravity of the fluid
    """

    volume_flowrate_gpm = convert_cubic_meters_per_second_to_gpm(
        volume_flowrate)
    return ((flow_coefficient / volume_flowrate_gpm) ** 2 / specific_gravity) ** -1


def calculate_volume_flowrate(mass_flowrate, density):
    """Calculates the volume flowrate from the mass flowrate and the density of the fluid.

    Keyword arguments:
    mass_flowrate   --  rate of flow of the fluid (kg/s)
    density         --  fluid density (kg/m^3)
    """

    return mass_flowrate / density


def calculate_pressure_after_valve(mass_flowrate, temperature, pressure, flow_coefficient):
    """Calculates the pressure in the fluid after traveling through the valve.

    Keyword arguments:
    mass_flowrate       --  rate of flow of the fluid (kg/s)
    temperature         --  temperature of the fluid (K)
    pressure            --  pressure of the fluid before the valve (Pa)
    flow_coefficient    --  flow coefficient of the valve
    """

    valve_density = CoolProp.CoolProp.PropsSI('D', 'T', temperature, 'P',
                                                   pressure, FLUID)
    valve_volume_flowrate = calculate_volume_flowrate(
        mass_flowrate, valve_density)
    valve_specific_gravity = calculate_specific_gravity(
        valve_density)
    valve_pressure_drop = calculate_flow_coefficient_pressure_drop(
        flow_coefficient, valve_volume_flowrate, valve_specific_gravity)
    return pressure - valve_pressure_drop


def calculate_total_pressure_drop(mass_flowrate, temperature, bottle_pressure):
    """Calculates the total pressure drop through the plumbing system.

    Keyword arguments:
    mass_flowrate   --  rate of flow of the fluid (kg/s)
    temperature     --  temperature of the fluid (K)
    bottle_pressure --  pressure in the bottle (Pa)
    """

    tank_valve_pressure = calculate_pressure_after_valve(
        mass_flowrate, temperature, bottle_pressure, TANK_VALVE_FLOW_COEFFICIENT)

    tee_1_pressure = tank_valve_pressure

    check_valve_pressure = calculate_pressure_after_valve(
        mass_flowrate, temperature, tee_1_pressure, CHECK_VALVE_FLOW_COEFFICIENT)

    tee_2_pressure = check_valve_pressure

    ball_valve_pressure = tee_2_pressure

    tee_3_pressure = ball_valve_pressure

    exit_pressure = tee_3_pressure

    print(bottle_pressure - exit_pressure)
    return bottle_pressure - exit_pressure


def main():
    # Surrounding atmospheric temperature (K)
    ATMOSPHERIC_TEMPERATURE = convert_f_to_k(ATMOSPHERIC_TEMPERATURE_F)

    # Mass of fluid in the bottle (kg)
    INITIAL_BOTTLE_MASS = convert_lb_to_kg(INITIAL_BOTTLE_WEIGHT_LBF)

    # Temperature inside the bottle (K)
    BOTTLE_TEMPERATURE = ATMOSPHERIC_TEMPERATURE

    # Gauge pressure inside the bottle (Pa)
    INITIAL_BOTTLE_GAUGE_PRESSURE = convert_psi_to_pa(
        INITIAL_BOTTLE_GAUGE_PRESSURE_PSI)

    # Absolute pressure inside the bottle (Pa)
    INITIAL_BOTTLE_ABSOLUTE_PRESSURE = convert_gauge_pressure_to_absolute(
        INITIAL_BOTTLE_GAUGE_PRESSURE, ATMOSPHERIC_PRESSURE)

    # Cross-sectional area of the pipe (m^2)
    PIPE_AREA = calculate_pipe_area(PIPE_INNER_DIAMETER)

    print("The bottle phase is", CoolProp.CoolProp.PhaseSI('P', INITIAL_BOTTLE_ABSOLUTE_PRESSURE, 'T', BOTTLE_TEMPERATURE, FLUID),
          "and the exit phase is", CoolProp.CoolProp.PhaseSI('P', ATMOSPHERIC_PRESSURE, 'T', ATMOSPHERIC_TEMPERATURE, FLUID))

    average_flow_pressure = (
        INITIAL_BOTTLE_ABSOLUTE_PRESSURE + ATMOSPHERIC_PRESSURE) / 2

    average_flow_temperature = (
        BOTTLE_TEMPERATURE + ATMOSPHERIC_TEMPERATURE) / 2

    average_flow_density = CoolProp.CoolProp.PropsSI(
        'D', 'T', average_flow_temperature, 'P', average_flow_pressure, FLUID)

    # Isentropic incompressible flow assumption
    isentropic_incompressible_flow_velocity = bernoulli_equation_exit_velocity(
        INITIAL_BOTTLE_ABSOLUTE_PRESSURE, 0, ATMOSPHERIC_PRESSURE, average_flow_density)

    isentropic_incompressible_mass_flowrate = calculate_mass_flowrate_from_velocity(
        isentropic_incompressible_flow_velocity, average_flow_density, PIPE_AREA)

    print("The isentropic incompressible flowrate is",
          round(convert_kg_per_s_to_lbm_per_s(isentropic_incompressible_mass_flowrate), 3), "lbm/s")

    # Incompressible flow assumption
    mass_flowrate_upper = isentropic_incompressible_mass_flowrate
    mass_flowrate_lower = 0.0000001

    pressure_drop_upper = calculate_total_pressure_drop(
        mass_flowrate_upper, BOTTLE_TEMPERATURE, INITIAL_BOTTLE_ABSOLUTE_PRESSURE)
    pressure_drop_lower = calculate_total_pressure_drop(
        mass_flowrate_lower, BOTTLE_TEMPERATURE, INITIAL_BOTTLE_ABSOLUTE_PRESSURE)

    real_pressure_drop = INITIAL_BOTTLE_ABSOLUTE_PRESSURE - ATMOSPHERIC_PRESSURE

    error_upper = pressure_drop_upper - real_pressure_drop
    error_lower = pressure_drop_lower - real_pressure_drop

    error_guess = 1

    # Bisection root-finding method
    while abs(error_guess) > (10 ** -3):
        mass_flowrate_guess = (mass_flowrate_lower + mass_flowrate_upper) / 2
        pressure_drop_guess = calculate_total_pressure_drop(
            mass_flowrate_guess, BOTTLE_TEMPERATURE, INITIAL_BOTTLE_ABSOLUTE_PRESSURE)
        error_guess = pressure_drop_guess - real_pressure_drop

        if error_lower * error_guess < 0:
            mass_flowrate_upper = mass_flowrate_guess
        elif error_lower * error_guess > 0:
            mass_flowrate_lower = mass_flowrate_guess

    incompressible_mass_flowrate = mass_flowrate_guess

    print("The incompressible flowrate is",
          round(incompressible_mass_flowrate, 3) / 0.45359237, "lbm/s")
    print("Pressure compare", round(pressure_drop_guess, 3) / (10 ** 3), "kPa to",
          round(INITIAL_BOTTLE_ABSOLUTE_PRESSURE - ATMOSPHERIC_PRESSURE, 3) / (10 ** 3), "kPa")


if __name__ == "__main__":
    main()

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
