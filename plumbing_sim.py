import CoolProp.CoolProp
import numpy
import math
import matplotlib.pyplot
import pandas
import scipy.optimize
import scipy.integrate
import scipy.stats
import seaborn

STEP_SIZE = 1e-0  # Spacing in seconds of recorded measurements
SEARCH_RESOLUTION = 1e-2  # Spacing in kg/s of mass flowrate estimation

FLUID = 'NITROUSOXIDE'

INITIAL_BOTTLE_WEIGHT_LBF = 10  # Weight of fluid in the bottle (lbf)
INITIAL_BOTTLE_GAUGE_PRESSURE_PSI = 650  # Pressure inside the bottle (psig)

ATMOSPHERIC_PRESSURE = 103081.67  # Surrounding atmospheric pressure (Pa)
ATMOSPHERIC_TEMPERATURE_F = 38  # Surrounding atmospheric temperature (F)

PIPE_INNER_DIAMETER = 0.00707  # Inner diamater of the pipe (m)

# Flow coefficients
# https://catalog.circlevalve.com/item/check-valves/200-series-check-valves-0-to-3000-psig/cs-259a-1pp
TANK_VALVE_FLOW_COEFFICIENT = 0.41
CHECK_VALVE_FLOW_COEFFICIENT = 1.6

# Loss coefficients
TEE_LOSS_COEFFICIENT = 0.9  # Line flow, threaded
BALL_VALVE_LOSS_COEFFICIENT = 0.05  # Fully open
EXIT_LOSS_COEFFICIENT = 1


def convert_f_to_k(temperature_f):
    """Converts temperature in degrees Fahrenheit into kelvins.

    Keyword arguments:
    temperature_f   -- temperature (f)
    """

    temperature = (temperature_f + 459.67) / 1.8
    return temperature


# Surrounding atmospheric temperature (K)
TEMPERATURE = convert_f_to_k(ATMOSPHERIC_TEMPERATURE_F)


def convert_lb_to_kg(mass_lb):
    """Converts mass in pounds into kilograms.

    Keyword arguments:
    mass_lb -- mass (lb)
    """

    mass = 0.45359237 * mass_lb
    return mass


# Mass of fluid in the bottle (kg)
INITIAL_BOTTLE_MASS = convert_lb_to_kg(INITIAL_BOTTLE_WEIGHT_LBF)


def convert_psi_to_pa(pressure_psi):
    """Converts pressure in pounds per square inch into pascals.

    Keyword arguments:
    pressure_psi    -- pressure (psi)
    """

    pressure = 6.894757e3 * pressure_psi
    return pressure


# Gauge pressure inside the bottle (Pa)
INITIAL_BOTTLE_GAUGE_PRESSURE = convert_psi_to_pa(
    INITIAL_BOTTLE_GAUGE_PRESSURE_PSI)


def convert_gauge_pressure_to_absolute(gauge_pressure):
    """Converts gauge pressure in pascals into absolute pressure in pascals.

    Keyword arguments:
    gauge_pressure  -- gauge pressure (Pa)
    """

    pressure = gauge_pressure + ATMOSPHERIC_PRESSURE
    return pressure


# Absolute pressure inside the bottle (Pa)
INITIAL_BOTTLE_PRESSURE = convert_gauge_pressure_to_absolute(
    INITIAL_BOTTLE_GAUGE_PRESSURE)

# Density of the fluid (kg/m^3)
INITIAL_FLUID_DENSITY = CoolProp.CoolProp.PropsSI(
    'D', 'P', INITIAL_BOTTLE_PRESSURE, 'T', TEMPERATURE, FLUID)

# Volume of the bottle (m^3)
BOTTLE_VOLUME = INITIAL_BOTTLE_MASS / INITIAL_FLUID_DENSITY


def calculate_pipe_area(inner_diameter):
    """Calculates the cross-sectional area in m^2 of the pipe.

    Keyword arguments:
    inner_diameter  -- inner diameter of the pipe (m)
    """

    pipe_area = math.pi * (inner_diameter / 2) ** 2
    return pipe_area


# Cross-sectional area of the pipe (m^2)
PIPE_AREA = calculate_pipe_area(PIPE_INNER_DIAMETER)


def bernoulli_equation_exit_velocity(pressure_in, velocity_in, pressure_out, density):
    """Calculates the exit flow velocity in m/s using the Bernoulli equation, assuming the fluid is incompressible and there are no friction losses.

    Keyword arguments:
    pressure_in -- entry pressure (Pa)
    velocity_in -- entry velocity (m/s)
    pressure_in -- exit pressure (Pa)
    density     -- density of the fluid (kg/m^3)
    """

    if not pressure_in > pressure_out:
        print(
            f"pressure_in, pressure_out ERROR IN bernoulli_equation_exit_velocity ({pressure_in}, {pressure_out})")
        return 0

    # The energy per unit mass of the fluid at the entry point
    energy_in = pressure_in + 1 / 2 * density * velocity_in ** 2

    # The velocity of the fluid at the exit point
    velocity_out = ((energy_in - pressure_out) * 2 / density) ** 0.5
    return velocity_out


def calculate_mass_flowrate_from_velocity(velocity, density):
    """Calculates the mass flowrate in kg/s of the fluid.

    Keyword arguments:
    velocity    -- flow velocity of the fluid (m/s)
    density     -- density of the fluid (kg/m^3)
    """

    mass_flowrate = PIPE_AREA * velocity * density
    return mass_flowrate


def convert_kg_per_s_to_lbm_per_s(mass_flowrate):
    """Converts mass flowrate in kg/s into lbm/s.

    Keyword arguments:
    mass_flowrate   -- mass flowrate of the fluid (kg/s)
    """

    mass_flowrate_lbm_per_s = mass_flowrate / 0.45359237
    return mass_flowrate_lbm_per_s


def calculate_volume_flowrate(mass_flowrate, density):
    """Calculates the volume flowrate from the mass flowrate and the density of the fluid.

    Keyword arguments:
    mass_flowrate   -- rate of flow of the fluid (kg/s)
    density         -- fluid density (kg/m^3)
    """

    volume_flowrate = mass_flowrate / density
    return volume_flowrate


def calculate_specific_gravity(density):
    """Calculates the specific gravity of a fluid.

    Keyword arguments:
    density -- fluid density (kg/m^3)
    """

    WATER_DENSITY = 1000  # (kg/m^3)

    specific_gravity = density / \
        WATER_DENSITY
    return specific_gravity


def convert_cubic_meters_per_second_to_gpm(volume_flowrate):
    """Converts volume flowrate in m^3/s into gal/min.

    Keyword arguments:
    volume_flowrate -- volume flowrate (m^3/s)
    """

    volume_flowrate_gpm = 1.585e4 * volume_flowrate
    return volume_flowrate_gpm


def calculate_flow_coefficient_pressure_drop(flow_coefficient, volume_flowrate, specific_gravity):
    """Calculates the pressure drop in Pa across a valve using the flow coefficient.

    Keyword arguments:
    flow_coefficient    -- flow coefficient of the valve
    volume_flowrate     -- rate of flow of the fluid (m^3/s)
    specific gravity    -- specific gravity of the fluid
    """

    if not volume_flowrate > 0:
        print(
            f"volume_flowrate ERROR IN calculate_flow_coefficient_pressure_drop ({volume_flowrate})")
        return 0

    volume_flowrate_gpm = convert_cubic_meters_per_second_to_gpm(
        volume_flowrate)

    pressure_drop_psi = (
        (flow_coefficient / volume_flowrate_gpm) ** 2 / specific_gravity) ** -1

    pressure_drop = convert_psi_to_pa(pressure_drop_psi)
    return pressure_drop


def calculate_valve_pressure_drop(mass_flowrate, pressure, flow_coefficient):
    """Calculates the pressure drop in pascals in the fluid while traveling through the valve.

    Keyword arguments:
    mass_flowrate       -- rate of flow of the fluid (kg/s)
    pressure            -- pressure of the fluid before the valve (Pa)
    flow_coefficient    -- flow coefficient of the valve
    """

    if not pressure > 0:
        print(f"pressure ERROR IN calculate_valve_pressure_drop ({pressure})")
        return 0

    try:
        density = CoolProp.CoolProp.PropsSI(
            'D', 'T', TEMPERATURE, 'P', pressure, FLUID)
    except:
        density = CoolProp.CoolProp.PropsSI(
            'D', 'T|liquid', TEMPERATURE, 'P', pressure, FLUID)

    valve_volume_flowrate = calculate_volume_flowrate(
        mass_flowrate, density)
    valve_specific_gravity = calculate_specific_gravity(
        density)

    valve_pressure_drop = calculate_flow_coefficient_pressure_drop(
        flow_coefficient, valve_volume_flowrate, valve_specific_gravity)
    return valve_pressure_drop


def calculate_velocity_from_mass_flowrate(mass_flowrate, density):
    """Calculates the velocity in m/s of the fluid.

    Keyword arguments:
    mass_flowrate   -- flowrate of the fluid (kg/s)
    density         -- density of the fluid (kg/m^3)
    """

    volume_flowrate = mass_flowrate / density

    velocity = volume_flowrate / PIPE_AREA
    return velocity


def calculate_pressure_drop_from_head_loss(mass_flowrate, pressure, head_loss_coefficient):
    """Calculates the pressure drop in pascals for a pipe component.

    Keyword arguments:
    mass_flowrate           -- rate of flow of the fluid (kg/s)
    pressure                -- pressure of the fluid before the valve (Pa)
    head_loss_coefficient   -- the loss coefficient for the pipe component
    """

    if not pressure > 0:
        print(
            f"pressure ERROR IN calculate_pressure_drop_from_head_loss ({pressure})")
        return 0

    try:
        density = CoolProp.CoolProp.PropsSI(
            'D', 'T', TEMPERATURE, 'P', pressure, FLUID)
    except:
        density = CoolProp.CoolProp.PropsSI(
            'D', 'T|liquid', TEMPERATURE, 'P', pressure, FLUID)

    flow_velocity = calculate_velocity_from_mass_flowrate(
        mass_flowrate, density)

    pressure_drop = head_loss_coefficient * 1 / 2 * density * flow_velocity ** 2
    return pressure_drop


def calculate_all_pressure_drops(mass_flowrate, bottle_pressure):
    """Calculates the pressure drop across each component.

    Keyword arguments:
    mass_flowrate   -- rate of flow of the fluid (kg/s)
    bottle_pressure -- pressure in the bottle (Pa)
    """

    defined = True

    if not mass_flowrate > 0:
        print(
            f"mass_flowrate ERROR IN calculate_all_pressure_drops ({mass_flowrate})")
        defined = False
        return {'defined': defined, 'total': 0, 'tank valve': 0, 'tee 1': 0, 'check valve': 0, 'tee 2': 0, 'ball valve': 0, 'tee 3': 0, 'exit': 0}

    if not bottle_pressure > 0:
        print(
            f"bottle_pressure ERROR IN calculate_all_pressure_drops ({bottle_pressure})")
        defined = False
        return {'defined': defined, 'total': 0, 'tank valve': 0, 'tee 1': 0, 'check valve': 0, 'tee 2': 0, 'ball valve': 0, 'tee 3': 0, 'exit': 0}

    tank_valve_pressure_drop = calculate_valve_pressure_drop(
        mass_flowrate, bottle_pressure, TANK_VALVE_FLOW_COEFFICIENT)
    tank_valve_pressure = bottle_pressure - tank_valve_pressure_drop

    if not tank_valve_pressure > 0:
        defined = False
        return {'defined': defined, 'total': 0, 'tank valve': tank_valve_pressure, 'tee 1': 0, 'check valve': 0, 'tee 2': 0, 'ball valve': 0, 'tee 3': 0, 'exit': 0}
    tee_1_pressure_drop = calculate_pressure_drop_from_head_loss(
        mass_flowrate, tank_valve_pressure, TEE_LOSS_COEFFICIENT)
    tee_1_pressure = tank_valve_pressure - tee_1_pressure_drop

    if not tee_1_pressure > 0:
        defined = False
        return {'defined': defined, 'total': 0, 'tank valve': tank_valve_pressure, 'tee 1': tee_1_pressure, 'check valve': 0, 'tee 2': 0, 'ball valve': 0, 'tee 3': 0, 'exit': 0}
    check_valve_pressure_drop = calculate_valve_pressure_drop(
        mass_flowrate, tee_1_pressure, CHECK_VALVE_FLOW_COEFFICIENT)
    check_valve_pressure = tee_1_pressure - check_valve_pressure_drop

    if not check_valve_pressure > 0:
        defined = False
        return {'defined': defined, 'total': 0, 'tank valve': tank_valve_pressure, 'tee 1': tee_1_pressure, 'check valve': check_valve_pressure, 'tee 2': 0, 'ball valve': 0, 'tee 3': 0, 'exit': 0}
    tee_2_pressure_drop = calculate_pressure_drop_from_head_loss(
        mass_flowrate, check_valve_pressure, TEE_LOSS_COEFFICIENT)
    tee_2_pressure = check_valve_pressure - tee_2_pressure_drop

    if not tee_2_pressure > 0:
        defined = False
        return {'defined': defined, 'total': 0, 'tank valve': tank_valve_pressure, 'tee 1': tee_1_pressure, 'check valve': check_valve_pressure, 'tee 2': tee_2_pressure, 'ball valve': 0, 'tee 3': 0, 'exit': 0}
    ball_valve_pressure_drop = calculate_pressure_drop_from_head_loss(
        mass_flowrate, tee_2_pressure, BALL_VALVE_LOSS_COEFFICIENT)
    ball_valve_pressure = tee_2_pressure - ball_valve_pressure_drop

    if not ball_valve_pressure > 0:
        defined = False
        return {'defined': defined, 'total': 0, 'tank valve': tank_valve_pressure, 'tee 1': tee_1_pressure, 'check valve': check_valve_pressure, 'tee 2': tee_2_pressure, 'ball valve': ball_valve_pressure, 'tee 3': 0, 'exit': 0}
    tee_3_pressure_drop = calculate_pressure_drop_from_head_loss(
        mass_flowrate, ball_valve_pressure, TEE_LOSS_COEFFICIENT)
    tee_3_pressure = ball_valve_pressure - tee_3_pressure_drop

    if not tee_3_pressure > 0:
        defined = False
        return {'defined': defined, 'total': 0, 'tank valve': tank_valve_pressure, 'tee 1': tee_1_pressure, 'check valve': check_valve_pressure, 'tee 2': tee_2_pressure, 'ball valve': ball_valve_pressure, 'tee 3': tee_3_pressure, 'exit': 0}
    exit_pressure_drop = calculate_pressure_drop_from_head_loss(
        mass_flowrate, tee_3_pressure, EXIT_LOSS_COEFFICIENT)
    exit_pressure = tee_3_pressure - exit_pressure_drop

    total_pressure_drop = bottle_pressure - exit_pressure

    return {'defined': True, 'total': total_pressure_drop, 'tank valve': tank_valve_pressure_drop, 'tee 1': tee_1_pressure_drop, 'check valve': check_valve_pressure_drop, 'tee 2': tee_2_pressure_drop, 'ball valve': ball_valve_pressure_drop, 'tee 3': tee_3_pressure_drop, 'exit': exit_pressure_drop}


'''
def calculate_total_pressure_drop_error(mass_flowrate, bottle_pressure, real_pressure_drop):
    """Calculates the error between real and guess in pressure drop through the plumbing system.

    Keyword arguments:
    mass_flowrate       -- rate of flow of the fluid (kg/s)
    bottle_pressure     -- pressure in the bottle (Pa)
    real_pressure_drop  -- the actual pressure drop through the plumbing system (Pa)
    """

    if not mass_flowrate > 0:
        print(
            f"mass_flowrate ERROR IN calculate_total_pressure_drop_error ({mass_flowrate})")
        return real_pressure_drop

    if not bottle_pressure > 0:
        print(
            f"bottle_pressure ERROR IN calculate_total_pressure_drop_error ({bottle_pressure})")
        return real_pressure_drop

    if not real_pressure_drop > 0:
        print(
            f"real_pressure_drop ERROR IN calculate_total_pressure_drop_error ({real_pressure_drop})")
        return real_pressure_drop

    estimated_pressure_drop = calculate_all_pressure_drops(
        mass_flowrate, bottle_pressure)['total']

    return estimated_pressure_drop - real_pressure_drop
'''


def calculate_incompressible_mass_flowrate(bottle_pressure):
    """Calculates the incompressible mass flowrate in m^3/s through the plumbing system.

    Keyword arguments:
    bottle_pressure    -- the absolute pressure inside the bottle (Pa)
    """

    if not bottle_pressure > 0:
        print(
            f"bottle_absolute_pressure ERROR IN calculate_incompressible_mass_flowrate ({bottle_pressure})")
        return 0

    real_pressure_drop = bottle_pressure - ATMOSPHERIC_PRESSURE

    average_flow_pressure = (
        bottle_pressure + ATMOSPHERIC_PRESSURE) / 2

    average_flow_density = CoolProp.CoolProp.PropsSI(
        'D', 'T', TEMPERATURE, 'P', average_flow_pressure, FLUID)

    # Isentropic incompressible flow assumption
    isentropic_incompressible_flow_velocity = bernoulli_equation_exit_velocity(
        bottle_pressure, 0, ATMOSPHERIC_PRESSURE, average_flow_density)

    isentropic_incompressible_mass_flowrate = calculate_mass_flowrate_from_velocity(
        isentropic_incompressible_flow_velocity, average_flow_density)

    mass_flowrate_upper = isentropic_incompressible_mass_flowrate
    mass_flowrate_lower = 0

    if not mass_flowrate_upper > mass_flowrate_lower:
        print(
            f"mass_flowrate_upper, mass_flowrate_lower ERROR IN calculate_incompressible_mass_flowrate ({mass_flowrate_upper}, {mass_flowrate_lower})")
        return 0

    estimate = None
    estimate_pressure_drop_error = None
    for i in numpy.arange(mass_flowrate_lower, mass_flowrate_upper + SEARCH_RESOLUTION, SEARCH_RESOLUTION):
        if i != 0:
            pressure_drops_i = calculate_all_pressure_drops(i, bottle_pressure)
            if pressure_drops_i['defined']:
                total_pressure_drop_error_i = abs(
                    pressure_drops_i['total'] - real_pressure_drop)
                if not estimate_pressure_drop_error or total_pressure_drop_error_i < estimate_pressure_drop_error:
                    estimate = i
                    estimate_pressure_drop_error = total_pressure_drop_error_i

    if estimate and estimate_pressure_drop_error:
        return estimate

    print("ERROR IN calculate_incompressible_mass_flowrate NO ESTIMATE FLOWRATE FOUND")
    return 0


def convert_pa_to_psi(pressure):
    """Converts pressure in pascals into pounds per square inch.

    Keyword arguments:
    pressure -- pressure (pa)
    """

    pressure_psi = pressure / 6.894757e3
    return pressure_psi


def convert_kg_to_lbm(mass):
    """Converts mass in kilograms into pounds.

    Keyword arguments:
    mass -- mass (kg)
    """

    mass_lb = mass / 0.45359237
    return mass_lb


def convert_absolute_pressure_to_gauge(pressure):
    """Converts absolute pressure in pascals into gauge pressure.

    Keyword arguments:
    pressure    -- pressure (Pa)
    """

    gauge_pressure = pressure - ATMOSPHERIC_PRESSURE
    return gauge_pressure


def calculate_bottle_pressure(bottle_mass):
    """Calculates the pressure inside the NOX bottle.

    Keyword arguments:
    bottle_mass -- the current mass of fluid in the bottle (kg)
    """

    if bottle_mass < 0:
        print(
            f"bottle_mass ERROR IN calculate_bottle_pressure ({bottle_mass})")
        return 0

    density = bottle_mass / BOTTLE_VOLUME

    bottle_pressure = CoolProp.CoolProp.PropsSI(
        'P', 'D', density, 'T', TEMPERATURE, FLUID)
    return bottle_pressure


def calculate_incompressible_mass_flowrate_from_mass(bottle_mass):
    """Calculates the incompressible mass flowrate through the plumbing system.

    Keyword arguments:
    bottle_mass -- the mass of fluid in the bottle (kg)
    """

    bottle_absolute_pressure = calculate_bottle_pressure(bottle_mass[0])
    return calculate_incompressible_mass_flowrate(bottle_absolute_pressure)


def main():
    print("The bottle phase is", CoolProp.CoolProp.PhaseSI('P', INITIAL_BOTTLE_PRESSURE, 'T', TEMPERATURE, FLUID),
          "and the exit phase is", CoolProp.CoolProp.PhaseSI('P', ATMOSPHERIC_PRESSURE, 'T', TEMPERATURE, FLUID))

    average_flow_pressure = (
        INITIAL_BOTTLE_PRESSURE + ATMOSPHERIC_PRESSURE) / 2

    average_flow_density = CoolProp.CoolProp.PropsSI(
        'D', 'T', TEMPERATURE, 'P', average_flow_pressure, FLUID)

    # Isentropic incompressible flow assumption
    isentropic_incompressible_flow_velocity = bernoulli_equation_exit_velocity(
        INITIAL_BOTTLE_PRESSURE, 0, ATMOSPHERIC_PRESSURE, average_flow_density)

    isentropic_incompressible_mass_flowrate = calculate_mass_flowrate_from_velocity(
        isentropic_incompressible_flow_velocity, average_flow_density)

    print("The isentropic incompressible flowrate is",
          round(convert_kg_per_s_to_lbm_per_s(isentropic_incompressible_mass_flowrate), 3), "lbm/s")

    # Incompressible flow assumption
    incompressible_mass_flowrate = calculate_incompressible_mass_flowrate(
        INITIAL_BOTTLE_PRESSURE)

    incompressible_pressure_drop = calculate_all_pressure_drops(
        incompressible_mass_flowrate, INITIAL_BOTTLE_PRESSURE)['total']

    print("The incompressible flowrate is", round(
        convert_kg_per_s_to_lbm_per_s(incompressible_mass_flowrate), 3), "lbm/s")
    print("Pressure compare", round(incompressible_pressure_drop, 3) / (10 ** 3), "kPa to",
          round(INITIAL_BOTTLE_PRESSURE - ATMOSPHERIC_PRESSURE, 3) / (10 ** 3), "kPa")

    final_time = INITIAL_BOTTLE_MASS / incompressible_mass_flowrate * 2
    solution = scipy.integrate.solve_ivp(fun=lambda t, y: -calculate_incompressible_mass_flowrate_from_mass(
        y), t_span=(0, final_time), y0=[INITIAL_BOTTLE_MASS], max_step=STEP_SIZE)

    bottle_mass_list = solution.y[0]

    print("IVP solved.")

    y_bottle_mass = []
    y_bottle_pressure = []
    y_mass_flowrate = []
    y_ball_valve_pressure_drop = []
    i = 0
    while i < len(solution.y[0]):
        i_mass_flowrate = calculate_incompressible_mass_flowrate_from_mass([
                                                                           bottle_mass_list[i]])

        y_mass_flowrate.append(convert_kg_to_lbm(i_mass_flowrate))

        y_bottle_mass.append(convert_kg_to_lbm(bottle_mass_list[i]))
        i_bottle_pressure = calculate_bottle_pressure(bottle_mass_list[i])
        y_bottle_pressure.append(convert_pa_to_psi(
            convert_absolute_pressure_to_gauge(i_bottle_pressure)))

        all_pressure_drops = calculate_all_pressure_drops(
            i_mass_flowrate, i_bottle_pressure)
        y_ball_valve_pressure_drop.append(
            convert_pa_to_psi(all_pressure_drops["ball valve"]))

        i = i + 1

    df = pandas.DataFrame(data={'Time (s)': solution.t, 'Bottle Fluid Mass (lbm)': y_bottle_mass, 'Mass Flowrate (lbm/s)': y_mass_flowrate,
                                'Bottle Gauge Pressure (psi)': y_bottle_pressure, 'Ball Valve Pressure Drop (psi)': y_ball_valve_pressure_drop})
    print(df)
    df.to_excel("output.xlsx")

    matplotlib.pyplot.figure()
    seaborn.lineplot(x='Time (s)', y='Bottle Fluid Mass (lbm)', data=df)

    matplotlib.pyplot.figure()
    seaborn.lineplot(x='Time (s)', y='Mass Flowrate (lbm/s)', data=df)

    matplotlib.pyplot.figure()
    seaborn.lineplot(x='Time (s)', y='Bottle Gauge Pressure (psi)', data=df)

    matplotlib.pyplot.figure()
    seaborn.lineplot(x="Time (s)", y='Ball Valve Pressure Drop (psi)', data=df)

    matplotlib.pyplot.show()


if __name__ == "__main__":
    main()
