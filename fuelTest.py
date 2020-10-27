# Helical Port Rocket simulator (Iterative Code)
# Written by: Thomas Kulin
# October 2, 2020

# Designed from Cole's matlab simulator and this study: https://arc-aiaa-org.proxy.queensu.ca/doi/pdf/10.2514/1.B36208

import numpy as np
import matplotlib.pyplot as plt
from sympy import nsolve, Symbol
from objects import CEA

# Variables for Hybrid Sim
finalD = 0.1143  # [m]    Maximum possible diameter of the motor after completed burn
initialD = 0.01524  # [m]    Initial diameter of combustion port, pre-burn
initial_OF = 5#1.1  # initial OF ratio
initial_Pc = 400  # initial chamber pressure [PSI]
r_helix = 0.00762*2  # [m] Radius of helix curvature
N_helix = 2#2.36  # number of helical port turns
Lp = 0.6  # [m]    Length of the combustion port
throatR = 0.035/2  # [m]    Radius of the nozzle throat
exitR = 0.041275  # [m]    Radius of the nozzle exit
lamda = 0.97  # Nozzle efficiency
Pa = 101325  # [Pa] Ambient pressure
R = 8314.41 / 29.19  # [J/kmol*K] Universal gas constant divided by MM of combustion products

mdot_ox = 1  # [kg/s] Oxidizer flow rate (experimentally measured)
MW_ox = 44.013  # [g/mol] oxidizer molecular weight
mu_ox = 2.7E-5  # absolute viscosity of N2O [(N*s)/m^2]. this value is for 20 C, but increases significantly with higher temperature
MW_fuel = 83.92  # [g/mol] ABS fuel molecular weight
rho_fuel = 975  # [kg/m^3]  Average density of ABS plastic
h_vap = 3  # [kJ/g] Heat of vaporization
T_vap = 600  # [K] vaporization temperature of the fuel. Estimate of the temperature at fuel surface for delta h calculation

# simulation parameters
dt = 0.05  # [s] Differential time step to be used for each iteration
maxIterations = 1000

# Caclulated initial variables for simulation
A_t = np.pi * pow(throatR, 2)  # [m^2]  Specify the nozzle throat area
A_e = np.pi * pow(exitR, 2)  # [m^2]  Calculate the nozzle exit area
epsilon = A_e / A_t  # Nozzle area expansion ratio


# Begin Simulation
r_0 = initialD / 2
r_final = finalD / 2
r_L = [0 for x in range(maxIterations)]
OF_ratio = [initial_OF for x in range(maxIterations)]
Pc = [initial_Pc for x in range(maxIterations)]
m_fuel = [0 for x in range(maxIterations)]
m_ox = [0 for x in range(maxIterations)]
G_total = [0 for x in range(maxIterations)]
rdot_helix_plot = [0 for x in range(maxIterations)]
mdot_fuel = [0 for x in range(maxIterations)]
mdot_total = [0 for x in range(maxIterations)]
thrust = [0 for x in range(maxIterations)]
time = [0 for x in range(maxIterations)]

test = [0 for x in range(maxIterations)]
r_L[0] = r_0
i = 0

# C = CEA(oxName='N2O', fuelName='ABS')
C = CEA(oxName='LO2', fuelName='ABS')

while r_L[i] < r_final and i < maxIterations-1:
    i=i+1  # step counter
    time[i] = time[i-1] + dt

    # Calculate thermodynamic properties in combustion chamber
    C.getOutput(Pc[i], OF_ratio[i], epsilon, True)
    Isp_Vac, C_star, T_flame, MW, gamma, Cp, mu, thermCond, prandtl = C.getChamberEquilibrium(Pc[i-1], OF_ratio[i-1], epsilon, 1)  # get ISP [s], C* [m/s], T [K], MW [g/mol], gamma, Cp [j/g*C], viscosity[Poise], thermal conductivity [W/m*K], Prandtl Number [-]
    Q_total = C.getReactionHeat(Pc[i-1], OF_ratio[i-1], epsilon, 1)
    deltaH_surf = abs(Q_total) - Cp*(T_vap - 298.15)/1000  # Convective enthalpy transfer per unit massflow from the flame zone to the fuel
    h_ratio = deltaH_surf / h_vap
    print("delta H surf / Hv  = ", h_ratio)
    test[i] = Cp

    """
     A_t : choked nozzle throat area [m^2]
     P_0 : combustion chamber pressure [kPa]
     gamma : ratio of specific heats (from RPA)
     R_g : gas constant for combustion products (from RPA) [J/(kg*K)]
     T_0 : combustion flame temperature (from RPA) [K] (T_Flame)
     rho_fuel : solid fuel grain material density [kg/(m^3)]
    """

    # regression calculations
    R_g = 8.314462 / MW * 1000  # Gas constant for combustion products [j/(kg*K)]

    G_ox = mdot_ox / (np.pi * r_L[i - 1] ** 2)  # [kg/(s*m^2)]

    thingy = 0.047/(prandtl**(2/3)*rho_fuel) * h_ratio**0.23 * (mu/Lp)**0.2
    rdot_straight = thingy * (G_ox**0.2 + 5/9*thingy * Lp/(2*r_L[i-1]))**4  # straight port regression rate [m/s]  (eqn 2.27 High Regression Rate Hybrid Rocket Fuel Grains with Helical Port Structures)

    mdot_fuel[i] = rdot_straight * rho_fuel * np.pi * r_L[i - 1] * 2 * Lp  # fuel mass flow rate [kg/s]

    r_L[i] = r_L[i - 1] + rdot_straight * dt  # mean of longitudinal port radius [cm]. may be wrong as per eqn. 35

    G_total_est = (mdot_fuel[i] + mdot_ox) / (np.pi * pow(r_L[i], 2))  # total mass flux estimate [kg/(m^2*s)]
    OF_est = mdot_ox / mdot_fuel[i]  # Oxidizer to fuel ratio estimate

    R_e = G_total_est * Lp / mu  # longitudinal Reynolds number

    C_f_straight = 0.074 / pow(R_e, 1/5)  # straight port skin friction coefficient

    P_helix = (Lp / N_helix) if N_helix is not 0 else 1e15  # [m] helix pitch distance
    R_c = r_helix * (1 + pow(P_helix / (2 * np.pi * r_helix), 2))  # helix radius of curvature [m]
    R_c_eff = R_c * np.sqrt(1 + np.pi / 2 * pow((2 * r_L[i] - 2 * r_0) / R_c, 2))  # effective helix radius of curvature [m] (accounts for change in port radius)

    A_C_f = 1 + 0.0075/C_f_straight * np.sqrt(r_L[i] / R_c_eff)  # skin friction amplification factor
    A_beta = (1 + (r_helix/r_L[i-1])/h_ratio)  # radial wall blowing amplification factor
    A_total = A_beta * A_C_f

    rdot_helix = A_C_f * A_beta * rdot_straight  # helical port regression rate [m/s]
    rdot_helix_plot[i] = rdot_helix

    m_fuel[i] = rdot_helix * rho_fuel * np.pi * r_L[i - 1] * 2 * Lp *dt  # fuel consumed in time step [kg]
    m_ox[i] = mdot_ox * dt  # oxidizer consumed in time step [kg]

    OF_ratio[i] = m_ox[i] / m_fuel[i]  # Oxidizer to fuel ratio
    # G_total[i] = (m_fuel[i]/dt + mdot_ox) / (np.pi * pow(r_L[i], 2))  # total mass flux estimate [kg/(m^2*s)]



    mdot_total[i] = (m_fuel[i] + m_ox[i]) / dt
    Pc_Pa = C_star * mdot_total[i] / (A_t)  # calculate combustion chamber pressure
    Pc[i] = Pc_Pa / 6894.76  # convert to PSI

    if i > 3:
        OF_ratio[i] = (OF_ratio[i]+OF_ratio[i-1]+OF_ratio[i-2])/3  # Take a moving average to negate OF oscillations
        Pc[i] = (Pc[i]+Pc[i-1]+Pc[i-2])/3
    Isp_Vac, C_star, T_flame, MW, gamma, Cp, mu, thermCond, prandtl = C.getChamberEquilibrium(Pc[i], OF_ratio[i], epsilon, 3)
    test[i] = gamma
    gamma = 1.17
    # Calcualte the exit Mach number, to do so use eqn 3.100 from SPAD (Humble)and solve it numerically
    x = Symbol('x')
    Me = float(nsolve(
        epsilon ** (2 * (gamma - 1) / (gamma + 1)) - (2 / (1 + gamma)) * x ** (-2 * (gamma - 1) / (gamma + 1)) - (
                    (gamma - 1) / 2) * x ** (2 * (1 - ((gamma - 1) / (gamma + 1)))), x, 3))  # initial guess 3

    # Calculate the exit pressure, use eqn 3.95 from SPAD (Humble), with the
    # previously determined chamber pressure as the stagnation pressure
    Pe = Pc_Pa * (1 + ((gamma - 1) / 2) * Me ** 2) ** -(gamma / (gamma - 1))  # [Pa]

    # Calculate the exit exhaust temperature to determine the exit velocity,
    # given by eqn 3.94 from SPAD (Humble)
    Te = ((1 + (((gamma - 1) / 2) * Me ** 2)) ** -1) * T_flame

    # Caclualte the exit velocity using eqn 3.112 from SPAD (Humble), using the
    # exit temperature from above
    Ve = np.sqrt(((2 * gamma * R_g * Te) / (gamma - 1)) * (1 - (Pe / Pc_Pa)) ** ((gamma - 1) / gamma))  # [m/s]

    # Now calculate the theoretical thrust of the motor using eqn 1.6 from SPAD
    thrust[i] = lamda * (mdot_total[i] * Ve + (Pe - Pa) * A_e)  # [N]

fig, ax = plt.subplots()
ax.plot(time, thrust, '.r', label='Thrust')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Thrust [N]')
ax.set_title('Thrust vs regression rate')
plt.grid()
ax2 = ax.twinx()
ax2.plot(time, rdot_helix_plot, '.g', label='rdot')
ax2.set_ylabel('Regression Rate [m/s]')
plt.grid()

fig, ax = plt.subplots()
ax.plot(time, OF_ratio, '.b', label='OF')
ax.set_xlabel('Time [s]')
ax.set_ylabel('OF Ratio [-]')
ax.set_title('OF vs. Time')
plt.grid()

fig, ax = plt.subplots()
ax.plot(time, Pc, '.', label='Pc')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Chamber Pressure [PSI]')
ax.set_title('Chamber Pressure vs. Time')
plt.grid()

fig, ax = plt.subplots()
ax.plot(time, mdot_total, '.r', label='mdot_fuel')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Total Mass Flow [kg/s]')
ax.set_title('Total Mass Flow Rate vs. Time')
plt.grid()

fig, ax = plt.subplots()
ax.plot(time, test, '.r', label='test')
ax.set_xlabel('Time [s]')
ax.set_ylabel('test')
ax.set_title('test')
plt.grid()

plt.show()
