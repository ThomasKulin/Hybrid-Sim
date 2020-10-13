# Classes and stuff
# Written by: Thomas Kulin
# October 2, 2020

from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import matlab.engine
import os


class Injector:
    def __init__(self, V=0.01824, M_0=13, T1_0=274.25, n_inj=20, d_inj=0.0015, Cd=0.8, P2=85.9e-3, species="N2O", dt=0.5):
        self.V = V  # Tank volume [m^3]
        self.M = M_0  # Tank fluid mass [kg]
        self.T1 = T1_0  # Tank fluid temperature [K]
        self.n_inj = n_inj  # Number of injector orifices
        self.d_inj = d_inj  # Injector orifice diameter [m]
        self.Cd = Cd  # Injector discharge coefficient
        self.P2 = P2  # downstream pressure [MPa]
        self.species = species  # species going through the injector (note, humans are not permitted to enter)
        self.dt = dt  # time step [s]
        self.initializeMatlab()

    def initializeMatlab(self):
        self.eng = matlab.engine.start_matlab()
        path = os.getcwd()
        self.eng.cd(path+'\\injectorSim')

    def simulate(self):
        print("INPUTS:", float(self.V), float(self.M), float(self.T1), float(self.n_inj), float(self.d_inj), float(self.Cd), float(self.P2), self.species, float(self.dt))
        data = self.eng.InjectorSim(float(self.V), float(self.M), float(self.T1), float(self.n_inj), float(self.d_inj), float(self.Cd), float(self.P2), self.species, float(self.dt), nargout=9)
        self.M = data[0]  # Tank fluid mass [kg]
        self.rho1 = data[1]  # Tank fluid density [kg/m^3]
        self.T1 = data[2]  # Tank Temperature [K]
        self.P1 = data[3]  # Tank pressure [MPa]
        self.X1 = data[4]  # Tank fluid quality (measurement of saturation)
        self.h1 = data[5]  # Tank specific enthalpy [kJ/kg]
        self.H1 = data[6]  # Tank total enthalpy [kJ]
        self.mdot = data[7]  # Tank mass flow rate [kg/s]
        self.state1 = data[8]  # Tank fluid State [-1=- Input, 0=Liq, ...1=Sat, 2=Gas]
        print(data)


class CEA:
    def __init__(self, oxName, fuelName, OFRatio, Pc, ExpansionRatio):
        self.C = CEA_Obj(oxName=oxName, fuelName=fuelName)
        output = self.C.get_full_cea_output(Pc=Pc, MR=OFRatio, eps=ExpansionRatio, short_output=1)

    def addABS(self):
        # Add ABS species to CEA
        card_str = """
        fuel ABS  C 3.85   H 4.85   N 0.43     wt%=100.00
        h,cal=14990    t(k)=298.15   rho=1
        """
        add_new_fuel('ABS', card_str)