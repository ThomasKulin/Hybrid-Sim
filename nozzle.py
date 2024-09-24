import math

# er = float(er)  # Expansion ratio
# ea = float(ea)  # exit area
# ta = float(ta)  # throat area
# mdot = float(mdot)  # Ideal propellant consumption in kg/s
tf = float(3500)  # Thrust force in N
k = float(1.17)  # gamma value=1.22
# pt = float(pt)  # critical pressure throat
# vt = float(vt)  # throat velocity
r = float(284.8376156)  # in J/(K*mol*K)
ct = float(3350)  # chamber temp in K
Pc_PSI = 400  # chamber pressure in PSI
p1 = float(Pc_PSI*0.006894757)  # chamber pressure in MPa
p1big = float(p1*1e6)  # ^^ in Pa
# iev = float(iev)  # ideal exit velocity
# v1 = float(v1)  # Specific volume v1 at nozzle entrance
# sv2 = float(sv2)  # Specific volume v2 at exit
# p2 = float(0.093217)  # presure at 0.75km in MPa
p2 = float(0.099566)    # pressure at 150m in MPa
# svt = float(svt)  # specific volume vt at throat
p = float(0.03377)  # pressure ratio

pt = p1*((2/(k+1))**(k/(k-1)))
print("critical pressure:", pt)

vt = math.sqrt(((2*k)/(k+1))*r*ct)
print("throat velocity:", vt)

iev = math.sqrt(((2*k/(k-1))*r*ct)*(1-(p**((k-1)/k))))
print("ideal exit velocity: ", iev)

mdot = tf/iev
print("Ideal consumption: ", mdot)

v1 = (r*ct)/p1big
print("Specific volume v1 at nozzle entrance:", v1)

svt = v1*(((k+1)/2)**(1/(k-1)))
print("specific volume vt at throat:", svt)

sv2 = v1*((p1/p2)**(1/k))
print("Specific volume v2 at exit:", sv2)

ta = (mdot*svt)/vt
print("throat area:", ta)

ea = (mdot*sv2)/iev
print("exit area:", ea)

er = ea/ta
print("expansion ratio:", er)