from objects import Injector
import matplotlib.pyplot as plt

ox_Species = "CO2"

# simulation parameters
dt = 0.0025  # [s] Differential time step to be used for each iteration
maxIterations = 2000

v_tank = 0.009  # [M^3] Oxidizer tank volume
m_tank = 0.17  # [kg] Mass of oxidizer in tank
t_tank = 273+14  # [K] Initial tank temperature
p_tank = 700  # [PSI] Initial tank pressure
n_inj = 26  # number of injector orifices
d_inj = 0.0015  # [m] Orifice port diameter
Cd = 0.8  # injector discharge coefficient
p_ambient = 14.6959  # Ambient pressure [PSI]

inj = Injector(r"C:\Users\Thomas\Documents\Projects\Propulsion 2020\Hybrid-Sim")

# Begin Simulation
m_ox = [0 for x in range(maxIterations)]
mdot_ox = [0 for x in range(maxIterations)]
temp = m_tank
m_tank = [temp for x in range(maxIterations)]
temp = t_tank
t_tank = [temp for x in range(maxIterations)]
temp = p_tank
p_tank = [temp for x in range(maxIterations)]
time = [0 for x in range(maxIterations)]
i = 0
while i < maxIterations-1:
    i=i+1  # step counter
    time[i] = time[i-1] + dt

    if (p_ambient > p_tank[i-1]-10):  # break out of sim if ambient pressure comes within 10 PSI of tank pressure. prevents errors at end of run
        break

    inj.initializeVariables(v_tank, m_tank[i - 1], t_tank[i - 1], p_tank[i-1] if p_tank[i-1] < 500 else 500, n_inj, d_inj, Cd, p_ambient / 145.038, ox_Species, dt)
    inj.simulate()
    mdot_ox[i] = inj.mdot
    m_tank[i] = inj.M
    t_tank[i] = inj.T1
    p_tank[i] = inj.P1 * 145.038  # convert MPA to PSI
    print("INJ PARAMS: ", "Flow:" + str(mdot_ox[i]) + "kg/s", "Tank:" + str(p_tank[i]) + 'PSI', "CO2 State:" + str(inj.state1))

    m_ox[i] = mdot_ox[i] * dt  # oxidizer consumed in time step [kg]


fig, ax = plt.subplots()
ax.plot(time, m_tank, '.b')
ax.set_xlabel('Time [s]')
ax.set_ylabel('CO2 Mass')
ax.set_title('CO2 Mass vs. Time')
plt.grid()
plt.legend()

fig, ax = plt.subplots()
ax.plot(time, [p_ambient for _ in range(len(time))], '.r', label='P_amb')
ax.plot(time, p_tank, '.b', label='P_tank')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Pressure [PSI]')
ax.set_title('Pressures vs. Time')
plt.grid()
plt.legend()

fig, ax = plt.subplots()
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass Flow [kg/s]')
ax.set_title('Flow Rate vs. Time')
plt.grid()
ax.plot(time, mdot_ox, '.b', label='mdot_ox')
plt.legend()

fig, ax = plt.subplots()
ax.plot(time, t_tank, '.r')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Temperature [K]')
ax.set_title('Tank Temperature vs. Time')
plt.grid()

plt.show()
print("done")
