# Helical Port Rocket simulator (Iterative Code)
# Written by: Thomas Kulin
# October 2, 2020

# Designed from Cole's matlab simulator and this study: https://arc-aiaa-org.proxy.queensu.ca/doi/pdf/10.2514/1.B36208

import numpy as np
import matplotlib.pyplot as plt
from sympy import nsolve, Symbol
from objects import CEA
from objects import Injector


# Variables for Hybrid Sim
v_tank = 0.01824  # [M^3] Oxidizer tank volume
m_tank = 13  # [kg] Mass of oxidizer in tank
t_tank = 298  # [K] Initial tank temperature
n_inj = 25  # number of injector orifices
d_inj = 0.0015  # [m] Orifice port diameter
Cd = 0.8  # injector discharge coefficient
ox_Species = "N2O"  # floaty boom stuff
initial_Pc = 50  # initial chamber pressure [PSI]
throatR = 0.035/2  # [m]    Radius of the nozzle throat
exitR = 0.041275  # [m]    Radius of the nozzle exit
lamda = 0.97  # Nozzle efficiency
Pa = 101325  # [Pa] Ambient pressure
R = 8314.41 / 29.19  # [J/kmol*K] Universal gas constant divided by MM of combustion products


# Oxidizer Parameters
# mdot_ox = 1  # [kg/s] Oxidizer flow rate (experimentally measured)
ox_Species = "N2O"
MW_ox = 44.013  # [g/mol] oxidizer molecular weight
mu_ox = 2.7E-5  # absolute viscosity of N2O [(N*s)/m^2]. this value is for 20 C, but increases significantly with higher temperature
initial_OF = 3#1.1  # initial OF ratio

# Fuel Parameters
finalD = 0.0254*5  # [m]    Maximum possible diameter of the motor after completed burn
initialD = 0.0254*2  # [m]    Initial diameter of combustion port, pre-burn
r_helix = 0.0254  # [m] Radius of helix curvature
N_helix = 1  # number of helical port turns
Lp = 0.4  # [m]    Length of the combustion port
MW_fuel = 83.92  # [g/mol] ABS fuel molecular weight
rho_fuel = 975  # [kg/m^3]  Average density of ABS plastic
h_vap = 3  # [kJ/g] Heat of vaporization
T_vap = 600  # [K] vaporization temperature of the fuel. Estimate of the temperature at fuel surface for delta h calculation

# simulation parameters
dt = 0.25  # [s] Differential time step to be used for each iteration
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
mdot_ox = [0 for x in range(maxIterations)]
mdot_fuel = [0 for x in range(maxIterations)]
mdot_total = [0 for x in range(maxIterations)]
temp = m_tank
m_tank = [temp for x in range(maxIterations)]
temp = t_tank
t_tank = [temp for x in range(maxIterations)]
p_tank = [750 for x in range(maxIterations)]
thrust = [0 for x in range(maxIterations)]
time = [0 for x in range(maxIterations)]

test = [0 for x in range(maxIterations)]
r_L[0] = r_0
i = 0

C = CEA(oxName=ox_Species, fuelName='ABS')
inj = Injector()

while r_L[i] < r_final and i < maxIterations-1 and Pc[i] > 40:
    i=i+1  # step counter
    time[i] = time[i-1] + dt

    if (Pc[i-1] > p_tank[i-1]-10):  # break out of sim if chamber pressure comes within 10 PSI of tank pressure. prevents errors at end of run
        break

    inj.initializeVariables(v_tank, m_tank[i-1], t_tank[i-1], n_inj, d_inj, Cd, Pc[i - 1] / 145.038, ox_Species, dt)
    inj.simulate()
    mdot_ox[i] = inj.mdot
    m_tank[i] = inj.M
    t_tank[i] = inj.T1
    p_tank[i] = inj.P1 * 145.038
    print("INJ PARAMS: ", "Flow:"+str(mdot_ox[i])+"kg/s", "Tank:"+str(p_tank[i])+'PSI', "Chamber:"+str(Pc[i-1])+"PSI")

    # Calculate thermodynamic properties in combustion chamber
    C.getOutput(Pc[i], OF_ratio[i], epsilon, True)
    Isp_Vac, C_star, T_flame, MW, gamma, Cp, mu, thermCond, prandtl = C.getChamberEquilibrium(Pc[i-1], OF_ratio[i-1], epsilon, 1)  # get ISP [s], C* [m/s], T [K], MW [g/mol], gamma, Cp [j/g*C], viscosity[Poise], thermal conductivity [W/m*K], Prandtl Number [-]
    Q_total = C.getReactionHeat(Pc[i-1], OF_ratio[i-1], epsilon, 1)
    deltaH_surf = abs(Q_total) - Cp*(T_vap - 298.15)/1000  # Convective enthalpy transfer per unit massflow from the flame zone to the fuel
    h_ratio = deltaH_surf / h_vap
    print("delta H surf / Hv  = ", h_ratio)


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

    G_ox = mdot_ox[i] / (np.pi * r_L[i - 1] ** 2)  # [kg/(s*m^2)]

    thingy = 0.047/(prandtl**(2/3)*rho_fuel) * h_ratio**0.23 * (mu/Lp)**0.2
    rdot_straight = thingy * (G_ox**0.2 + 5/9*thingy * Lp/(2*r_L[i-1]))**4  # straight port regression rate [m/s]  (eqn 2.27 High Regression Rate Hybrid Rocket Fuel Grains with Helical Port Structures)

    mdot_fuel[i] = rdot_straight * rho_fuel * np.pi * r_L[i - 1] * 2 * Lp  # fuel mass flow rate [kg/s]

    r_L[i] = r_L[i - 1] + rdot_straight * dt  # mean of longitudinal port radius [cm]. may be wrong as per eqn. 35

    G_total_est = (mdot_fuel[i] + mdot_ox[i]) / (np.pi * pow(r_L[i], 2))  # total mass flux estimate [kg/(m^2*s)]
    OF_est = mdot_ox[i] / mdot_fuel[i]  # Oxidizer to fuel ratio estimate

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
    mdot_fuel[i] = m_fuel[i]/dt
    m_ox[i] = mdot_ox[i] * dt  # oxidizer consumed in time step [kg]
    mdot_total[i] = (m_fuel[i] + m_ox[i]) / dt

    OF_ratio[i] = m_ox[i] / m_fuel[i]  # Oxidizer to fuel ratio

    Pc_Pa = C_star * mdot_total[i] / (A_t)  # calculate combustion chamber pressure
    Pc[i] = Pc_Pa / 6894.76  # convert to PSI

    if i > 3:
        OF_ratio[i] = (OF_ratio[i]+OF_ratio[i-1]+OF_ratio[i-2])/3  # Take a moving average to negate OF oscillations
        Pc[i] = (Pc[i]+Pc[i-1]+Pc[i-2])/3
    else:
        OF_ratio[i] = (OF_ratio[i]+OF_ratio[i-1])/2
        Pc[i] = (Pc[i]+Pc[i-1])/3
    Isp_Vac, C_star, T_flame, MW, gamma, Cp, mu, thermCond, prandtl = C.getChamberEquilibrium(Pc[i], OF_ratio[i], epsilon, 3)
    test[i] = gamma
    # gamma = 1.17
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
ax.set_ylim(0, 6000)
plt.grid()
ax2 = ax.twinx()
ax2.plot(time, rdot_helix_plot, '.g', label='rdot')
ax2.set_ylabel('Regression Rate [m/s]')
ax2.set_ylim(0, 0.08)
plt.grid()

fig, ax = plt.subplots()
ax.plot(time, OF_ratio, '.b', label='OF')
ax.set_xlabel('Time [s]')
ax.set_ylabel('OF Ratio [-]')
ax.set_title('OF vs. Time')
plt.grid()

fig, ax = plt.subplots()
ax.plot(time, Pc, '.r', label='Pc')
ax.plot(time, p_tank, '.b', label='P_tank')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Pressure [PSI]')
ax.set_title('Pressures vs. Time')
plt.grid()
plt.legend()

fig, ax = plt.subplots()
ax.plot(time, mdot_fuel, '.r', label='mdot_fuel')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass Flow [kg/s]')
ax.set_title('Flow Rate vs. Time')
plt.grid()
ax.plot(time, mdot_ox, '.b', label='mdot_ox')
ax.set_ylim(0, 2)
plt.legend()

fig, ax = plt.subplots()
ax.plot(time, test, '.r', label='test')
ax.set_xlabel('Time [s]')
ax.set_ylabel('test')
ax.set_title('test')
plt.grid()

plt.show()
