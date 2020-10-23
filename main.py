# Helical Port Rocket simulator (Iterative Code)
# Written by: Thomas Kulin
# October 2, 2020

# Designed from Cole's matlab simulator and this study: https://arc-aiaa-org.proxy.queensu.ca/doi/pdf/10.2514/1.B36208

import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant
from sympy import nsolve, Symbol

# Variables for Hybrid Sim
finalD = 0.1143     #[m]    Maximum possible diameter of the motor after completed burn
initialD = 0.0635  #[m]    Initial diameter of combustion port, pre-burn
r_helix = 0.02  # [m] Radius of helix curvature
N_helix = 3  # number of helical port turns
Lp = 0.60         # [m]    Length of the combustion port
mdot_ox = 1         # [kg/s] Oxidizer flow rate (experimentally measured)
MW_ox = 44.013  # [g/mol] oxidizer molecular weight
throatR = 0.018415    # [m]    Radius of the nozzle throat
exitR = 0.041275     # [m]    Radius of the nozzle exit
initial_OF = 3  # initial OF ratio
initial_Pc =400  # initial chamber pressure [PSI]

# Caclulated initial variables for simulation
A_t = np.pi*pow(throatR, 2)  # [m^2]  Specify the nozzle throat area
A_e = np.pi*pow(exitR, 2)   # [m^2]  Calculate the nozzle exit area
epsilon = A_e/A_t    # Nozzle area expansion ratio

rho_fuel = 975           # [kg/m^3]  Average density of ABS plastic
mu_ox = 1.47E-5  # absolute viscosity of N2O [(N*s)/m^2]. this value is for 20 C, but increases significantly with higher temperature
lamda = 0.97        # Nozzle efficiency
Pa = 101325         # [Pa] Ambient pressure
R = 8314.41/29.19   # [J/kmol*K] Universal gas constant divided by MM of combustion products

# simulation parameters
dt = 0.01            # [s] Differential time step to be used for each iteration
maxIterations = 10

# Begin Simulation
r_0 = initialD/2
r_final = finalD/2
r_L = [0 for x in range(maxIterations)]
OF_ratio = [initial_OF for x in range(maxIterations)]
Pc = [initial_Pc for x in range(maxIterations)]
m_fuel = [0 for x in range(maxIterations)]
m_ox = [0 for x in range(maxIterations)]
mdot_fuel = [0 for x in range(maxIterations)]
mdot_total = [0 for x in range(maxIterations)]
thrust = [0 for x in range(maxIterations)]
time = [0 for x in range(maxIterations)]
r_L[0] = r_0
i=0

# Add ABS species to CEA
card_str = """
       fuel ABS  C 3.85   H 4.85   N 0.43     wt%=100.00
       h,cal=14990    t(k)=298.15   rho=1
       """
add_new_fuel('ABS', card_str)

C = CEA_Obj(oxName='N2O', fuelName='ABS')


while r_L[i] < r_final and i < maxIterations-1:
    i=i+1  # step counter
    time[i] = time[i-1] + dt


    output = C.get_full_cea_output(Pc=Pc[i-1], MR=OF_ratio[i-1], eps=epsilon, short_output=0, output='siunits')
    #print(output)

    Isp_Vac, C_star, T_flame, MW, gamma = C.get_IvacCstrTc_ChmMwGam(Pc=Pc[i-1], MR=OF_ratio[i-1], eps=epsilon)  # get ISP [s], C* [ft/s], T [K], MW [g/mol]
    C_star = C_star/3.28084  # convert to [m/s]


    """
     A_t : choked nozzle throat area [m^2]
     P_0 : combustion chamber pressure [kPa] 
     gamma : ratio of specific heats (from RPA)
     R_g : gas constant for combustion products (from RPA) [J/(kg*K)]
     T_0 : combustion flame temperature (from RPA) [K] (T_Flame)
     rho_fuel : solid fuel grain material density [kg/(m^3)]
    """
    Pc_Pa = Pc[i-1] * 6894.76  # convert PSI to Pa

    # regression calculations
    R_g = 8.314462 / MW * 1000  # Gas constant for combustion products [j/(kg*K)]
    #mdot_total[i] = A_t * Pc_Pa/1000 * np.sqrt(gamma / (R_g * T_flame) * pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1))) / 1000  # nozzle exit mass flow rate [kg/s]

    alpha = 0.000007  # regression rate burn model coefficient [cm∕s]
    n = 0.8
    m = -0.2
    G_ox_est = mdot_ox/(np.pi * r_L[i-1]**2)  # [kg/(s*m^2)]
    rdot_straight_est = alpha * pow(G_ox_est, n) * pow(Lp, m)  # straight port regression rate starting point [m/s]

    mdot_fuel_est = rdot_straight_est * rho_fuel * np.pi * r_L[i-1] * 2 * Lp  # fuel mass flow rate [kg/s]


    #rdot_mean = mdot_fuel / (2 * np.pi * rho_fuel * r_L[i-1] * Lp)  # mean longitudinal regression rate [cm/s]
    r_L[i] = r_L[i-1] + rdot_straight_est * dt  # mean of longitudinal port radius [cm]. may be wrong as per eqn. 35
    #r_L[i] = np.sqrt(pow(r_0, 2) + 1/(np.pi*rho_fuel*Lp) * mdot_fuel[i] * time[i])  # instantaneous mean of longitudinal port radius [cm]


    G_total_est = (mdot_fuel_est + mdot_ox) / (np.pi * pow(r_L[i], 2))  # total mass flux estimate [kg/(m^2*s)]
    OF_est = mdot_ox / mdot_fuel_est  # Oxidizer to fuel ratio estimate

    R_e = G_total_est * Lp / mu_ox  # longitudinal Reynolds number

    C_f_straight = 0.074/pow(R_e, 1/5)  # straight port skin friction coefficient

    P_helix = Lp/N_helix  # [m] helix pitch distance
    R_c = r_helix * (1 + pow(P_helix/(2*np.pi*r_helix), 2))  # helix radius of curvature [m]
    R_c_eff = R_c * np.sqrt(1 + np.pi/2 * pow((2*r_L[i]-2*r_0)/R_c, 2))  # effective helix radius of curvature [m] (accounts for change in port radius)

    A_C_f = 1 + 0.0075/C_f_straight * np.sqrt(r_L[i]/R_c_eff)  # radial wall blowing amplification factor
    A_beta = pow(1 + 2 * pow(2*np.pi*N_helix, 2) * pow(OF_est, 2) * MW / MW_ox * r_helix / r_L[i], 0.77)  # radial wall blowing amplification factor

    rdot_straight = alpha * pow(G_total_est, n) * pow(Lp, m)  # straight port regression rate [m/s]
    rdot_helix = A_C_f * A_beta * rdot_straight  # helical port regression rate [m/s]

    #m_fuel[i] = 2*np.pi*Lp*rho_fuel*r_L[i]*rdot_helix*dt  # fuel consumed in time step [kg]
    m_fuel[i] = rdot_helix*dt
    m_ox[i] = mdot_ox * dt  # oxidizer consumed in time step [kg]

    OF_ratio[i] = m_ox[i]/m_fuel[i]  # Oxidizer to fuel ratio



    mdot_total[i] = (m_fuel[i] + m_ox[i])/dt
    Pc_Pa = C_star*mdot_total[i]/(A_t*0.0001)  # calculate combustion chamber pressure
    Pc[i] = Pc_Pa / 6894.76  # convert to PSI

    # Calcualte the exit Mach number, to do so use eqn 3.100 from SPAD (Humble)and solve it numerically
    x = Symbol('x')
    Me = float(nsolve(epsilon**(2*(gamma-1)/(gamma+1)) - (2/(1+gamma)) * x**(-2*(gamma-1)/(gamma+1)) - ((gamma-1)/2) * x**(2*(1-((gamma-1)/(gamma+1)))), x, 3))  # initial guess 3

    # Calculate the exit pressure, use eqn 3.95 from SPAD (Humble), with the
    # previously determined chamber pressure as the stagnation pressure
    Pe = Pc_Pa*(1+((gamma-1)/2)*Me**2)**-(gamma/(gamma-1))  # [Pa]
    
    # Calculate the exit exhaust temperature to determine the exit velocity,
    # given by eqn 3.94 from SPAD (Humble)
    Te = ((1 + (((gamma-1)/2)*Me**2))**-1)*T_flame
    
    # Caclualte the exit velocity using eqn 3.112 from SPAD (Humble), using the
    # exit temperature from above
    Ve = np.sqrt(((2*gamma*R_g*Te)/(gamma-1))*(1-(Pe/Pc_Pa))**((gamma-1)/gamma))  # [m/s]
    
    # Now calculate the theoretical thrust of the motor using eqn 1.6 from SPAD
    thrust[i] = lamda*(mdot_total[i]*Ve + (Pe-Pa)*A_e)  # [N]



plt.plot(thrust)
plt.show()
