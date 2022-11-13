import CoolProp.CoolProp
import math
import scipy.optimize
import scipy.integrate
import pandas
import seaborn
import matplotlib.pyplot

FLUID = 'NITROUSOXIDE'

INITIAL_BOTTLE_WEIGHT_LBF = 10  # Weight of fluid in the bottle (lbf)
INITIAL_BOTTLE_GAUGE_PRESSURE_PSI = 950  # Pressure inside the bottle (psig)

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


def convert_kg_to_lbm(mass_kg):
    """Converts mass in kilograms into pounds.

    Keyword arguments:
    mass_kg --  mass (kg)
    """

    return mass_kg / 0.45359237


def convert_psi_to_pa(pressure_psi):
    """Converts pressure in pounds per square inch into pascals.

    Keyword arguments:
    pressure_psi    --  pressure (psi)
    """

    return 6.894757e3 * pressure_psi


def convert_gauge_pressure_to_absolute(pressure):
    """Converts gauge pressure in pascals into absolute pressure.

    Keyword arguments:
    pressure                --  pressure (Pa)
    atmospheric_pressure    --  atmospheric pressure (Pa)
    """

    return pressure + ATMOSPHERIC_PRESSURE


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
    if not pressure_in > pressure_out:
        return 0

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

    if not volume_flowrate > 0:
        return 0
    volume_flowrate_gpm = convert_cubic_meters_per_second_to_gpm(
        volume_flowrate)
    return convert_psi_to_pa(((flow_coefficient / volume_flowrate_gpm) ** 2 / specific_gravity) ** -1)


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

    if not pressure > 0:
        return 0
    valve_density = CoolProp.CoolProp.PropsSI('D', 'T|gas', temperature, 'P', pressure, FLUID)
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

    return bottle_pressure - exit_pressure


def calculate_total_pressure_drop_error(mass_flowrate, temperature, bottle_pressure, real_pressure_drop):
    """Calculates the error between real and guess in pressure drop through the plumbing system.

    Keyword arguments:
    mass_flowrate   --  rate of flow of the fluid (kg/s)
    temperature     --  temperature of the fluid (K)
    bottle_pressure --  pressure in the bottle (Pa)
    real_pressure_drop  --  the actual pressure drop through the plumbing system (Pa)
    """

    return calculate_total_pressure_drop(mass_flowrate, temperature, bottle_pressure) - real_pressure_drop


def calculate_incompressible_mass_flowrate(bottle_absolute_pressure, temperature, pipe_area):
    """Calculates the incompressible mass flowrate through the plumbing system.

    Keyword arguments:
    bottle_absolute_pressure    --  the absolute pressure inside the bottle (Pa)
    temperature          --  the temperature in the bottle (K)
    pipe_area                   --  the area of the pipe (m^2)
    """

    real_pressure_drop = bottle_absolute_pressure - ATMOSPHERIC_PRESSURE

    average_flow_pressure = (
        bottle_absolute_pressure + ATMOSPHERIC_PRESSURE) / 2

    average_flow_density = CoolProp.CoolProp.PropsSI(
        'D', 'T', temperature, 'P', average_flow_pressure, FLUID)

    # Isentropic incompressible flow assumption
    isentropic_incompressible_flow_velocity = bernoulli_equation_exit_velocity(
        bottle_absolute_pressure, 0, ATMOSPHERIC_PRESSURE, average_flow_density)

    isentropic_incompressible_mass_flowrate = calculate_mass_flowrate_from_velocity(
        isentropic_incompressible_flow_velocity, average_flow_density, pipe_area)

    mass_flowrate_upper = isentropic_incompressible_mass_flowrate
    mass_flowrate_lower = 0

    if not mass_flowrate_upper > mass_flowrate_lower:
        return 0

    return scipy.optimize.bisect(
        calculate_total_pressure_drop_error, mass_flowrate_lower, mass_flowrate_upper, args=(temperature, bottle_absolute_pressure, real_pressure_drop))


def calculate_bottle_pressure(bottle_volume, bottle_mass, temperature):
    """Calculates the pressure inside the NOX bottle.

    Keyword arguments:
    bottle_volume   --  the internal volume of the bottle (m^3)
    bottle_mass     --  the current mass of fluid in the bottle (kg)
    temperature     --  the temperature in the bottle (K)
    """

    if bottle_mass < 0:
        return 0
    density = bottle_mass / bottle_volume
    return CoolProp.CoolProp.PropsSI(
        'P', 'D', density, 'T', temperature, FLUID)


def calculate_incompressible_mass_flowrate_from_mass(bottle_mass, bottle_volume, temperature, pipe_area):
    """Calculates the incompressible mass flowrate through the plumbing system.

    Keyword arguments:
    bottle_mass                 --  the mass of fluid in the bottle (kg)
    temperature                 --  the temperature in the bottle (K)
    pipe_area                   --  the area of the pipe (m^2)
    """

    bottle_absolute_pressure = calculate_bottle_pressure(bottle_volume,
                                                         bottle_mass[0], temperature)
    return calculate_incompressible_mass_flowrate(bottle_absolute_pressure, temperature, pipe_area)


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
        INITIAL_BOTTLE_GAUGE_PRESSURE)

    BOTTLE_VOLUME = INITIAL_BOTTLE_MASS / CoolProp.CoolProp.PropsSI(
        'D', 'P', INITIAL_BOTTLE_ABSOLUTE_PRESSURE, 'T', BOTTLE_TEMPERATURE, FLUID)

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
    incompressible_mass_flowrate = calculate_incompressible_mass_flowrate(
        INITIAL_BOTTLE_ABSOLUTE_PRESSURE, BOTTLE_TEMPERATURE, PIPE_AREA)

    incompressible_pressure_drop = calculate_total_pressure_drop(
        incompressible_mass_flowrate, BOTTLE_TEMPERATURE, INITIAL_BOTTLE_ABSOLUTE_PRESSURE)

    print("The incompressible flowrate is", round(
        convert_kg_per_s_to_lbm_per_s(incompressible_mass_flowrate), 3), "lbm/s")
    print("Pressure compare", round(incompressible_pressure_drop, 3) / (10 ** 3), "kPa to",
          round(INITIAL_BOTTLE_ABSOLUTE_PRESSURE - ATMOSPHERIC_PRESSURE, 3) / (10 ** 3), "kPa")

    final_time = INITIAL_BOTTLE_MASS / incompressible_mass_flowrate * 1.5
    solution = scipy.integrate.solve_ivp(fun=lambda t, y: -calculate_incompressible_mass_flowrate_from_mass(
        y, BOTTLE_VOLUME, BOTTLE_TEMPERATURE, PIPE_AREA), t_span=(0, final_time), y0=[INITIAL_BOTTLE_MASS])

    print(solution.t)
    print(solution.y[0])

    y_lbm = []
    for mass in solution.y[0]:
        y_lbm.append(convert_kg_to_lbm(mass))

    seaborn.set()
    seaborn.lineplot(x=solution.t, y=y_lbm)
    matplotlib.pyplot.show()


if __name__ == "__main__":
    main()