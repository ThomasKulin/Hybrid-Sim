# Classes and stuff
# Written by: Thomas Kulin
# October 2, 2020

from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import matlab.engine
import os
import numpy as np


class Injector:
    def __init__(self, path=None):
        self.V = 0  # Tank volume [m^3]
        self.M = 0  # Tank fluid mass [kg]
        self.T1 = 0 # Tank fluid temperature [K]
        self.P1 = 750 / 145.038  # initial tank pressure [mPa]
        self.n_inj = 1  # Number of injector orifices
        self.d_inj = 0.001  # Injector orifice diameter [m]
        self.Cd = 0  # Injector discharge coefficient
        self.P2 = 0  # downstream pressure [MPa]
        self.species = 0  # species going through the injector (note, humans are not permitted to enter)
        self.dt = 0  # time step [s]
        self.rho1 = 0  # Fluid density [kg/m^3]
        self.X1 = 0  # Fluid quality (measurement of saturation)
        self.h1 = 0  # Tank specific enthalpy [kJ/kg]
        self.H1 = 0  # Tank total enthalpy [kJ/kg]
        self.state1 = 0   # Tank fluid State [-1=- Input, 0=Liq, ...1=Sat, 2=Gas]
        self.initializeMatlab(path)

    def initializeVariables(self, V=0.01824, M_0=13, T1_0=274.25, P1_0=5.171, n_inj=20, d_inj=0.0015, Cd=0.8, P2=85.9e-3, species="N2O", dt=0.5):
        self.V = V  # Tank volume [m^3]
        self.M = M_0  # Tank fluid mass [kg]
        self.T1 = T1_0  # Tank fluid temperature [K]
        self.P1 = P1_0   # initial tank pressure [mPa]
        self.n_inj = n_inj  # Number of injector orifices
        self.d_inj = d_inj  # Injector orifice diameter [m]
        self.Cd = Cd  # Injector discharge coefficient
        self.P2 = P2  # downstream pressure [MPa]
        self.species = species  # species going through the injector (note, humans are not permitted to enter)
        self.dt = dt  # time step [s]

    def initializeMatlab(self, path=None):
        self.eng = matlab.engine.start_matlab()
        if not path:
            path = os.getcwd()
            self.eng.cd(path+'\\injectorSim')
        else:
            self.eng.cd(path + '\\injectorSim')


    def setChamberPressure(self, Pc):
        self.P2 = Pc

    def simulate(self):
        print("INPUTS:", float(self.V), float(self.M), float(self.T1), float(self.n_inj), float(self.d_inj), float(self.Cd), float(self.P2), self.species, float(self.dt))
        if (self.P1 > self.P2):
            data = self.eng.InjectorSim(float(self.V), float(self.M), float(self.T1), float(self.n_inj), float(self.d_inj), float(self.Cd), float(self.P2), self.species, float(self.dt), nargout=9)
        else:
            print("ERROR: Chamber pressure greater than tank pressure.")
            data = [self.M, self.rho1, self.T1, self.P1, self.X1, self.h1, self.H1, 0, self.state1]
            raise Exception("encountered situation that makes no sense")
        self.M = data[0]  # Tank fluid mass [kg]
        self.rho1 = data[1]  # Tank fluid density [kg/m^3]
        self.T1 = data[2]  # Tank Temperature [K]
        self.P1 = data[3]  # Tank pressure [MPa]
        self.X1 = data[4]  # Tank fluid quality (measurement of saturation)
        self.h1 = data[5]  # Tank specific enthalpy [kJ/kg]
        self.H1 = data[6]  # Tank total enthalpy [kJ]
        self.mdot = data[7]  # Tank mass flow rate [kg/s]
        self.state1 = data[8]  # Tank fluid State [-1=- Input, 0=Liq, ...1=Sat, 2=Gas]
        #print(data)


class Fuel:
    def __init__(self, r_L, mdot_ox):

        self.r_L = r_L
        initial_OF = 5  # 1.1  # initial OF ratio
        initial_Pc = 400  # initial chamber pressure [PSI]
        r_helix = 0.00762  # [m] Radius of helix curvature
        N_helix = 2.36  # number of helical port turns
        Lp = 0.3598  # [m]    Length of the combustion port
        throatR = 0.018415  # [m]    Radius of the nozzle throat
        exitR = 0.041275  # [m]    Radius of the nozzle exit
        lamda = 0.97  # Nozzle efficiency
        Pa = 101325  # [Pa] Ambient pressure
        R = 8314.41 / 29.19  # [J/kmol*K] Universal gas constant divided by MM of combustion products

        self.mdot_ox = mdot_ox  # [kg/s] Oxidizer flow rate (experimentally measured)
        MW_ox = 44.013  # [g/mol] oxidizer molecular weight
        mu_ox = 2.7E-5  # absolute viscosity of N2O [(N*s)/m^2]. this value is for 20 C, but increases significantly with higher temperature
        MW_fuel = 83.92  # [g/mol] ABS fuel molecular weight
        rho_fuel = 975  # [kg/m^3]  Average density of ABS plastic
        h_vap = 3  # [kJ/g] Heat of vaporization
        T_vap = 600  # [K] vaporization temperature of the fuel. Estimate of the temperature at fuel surface for delta h calculation

        # simulation parameters
        dt = 0.01  # [s] Differential time step to be used for each iteration
        maxIterations = 1000

        # Caclulated initial variables for simulation
        A_t = np.pi * pow(throatR, 2)  # [m^2]  Specify the nozzle throat area
        A_e = np.pi * pow(exitR, 2)  # [m^2]  Calculate the nozzle exit area
        epsilon = A_e / A_t  # Nozzle area expansion ratio


class CEA:
    def __init__(self, oxName, fuelName):
        self.oxName = oxName
        self.fuelName = fuelName
        # dictionary containing the standard heat of formation for all species in [kJ/mol]
        # you can find this data in thermo.inp in the rocketCEA library folder
        self.heatDict = {
            "ABS": 62.719,
            "Acrylic": 382.0,
            "N2O": 82.05,
            "LO2": 0,
            "*CO": -110.53,
            "*CO2": -393.5,
            "*H": 217.998,
            "HCO": 42.397,
            "HO2": 12.02,
            "*H2": 0,
            "C(gr)": 0,
            "H2O": -241.826,
            "*N": 472.680,
            "*NO": 91.271,
            "*N2": 0,
            "*O": 249.175,
            "*OH": 37.278,
            "*O2": 0,
            "COOH": -213,
            "H2O2": -135.88,
            "O3": 141.8,
            "CH3": 146.658,
            "CH4": -74.6,
            "C2H2,acetylene": 228.2,
            "CH2CO,ketene": -49.576,
            "C2H4": 52.5,
            "HCN": 133.082,
            "HNC": 194.378,
            "NH3": -54.752,
            "CH3CN": 31.38,
            "C4H2,butadiyne": 450,
            "C2N2": 283.209,
            "*CN": 438.683,
            "HNCO": -118.056,
            "C3H3,2-propynl": 331.8,
            "HNO": 102.032,
            "NO2": 34.193,
            "C2H6": -103.819,
            "HCHO,formaldehy": -108.58,
            "C3H6,propylene": 20.0
        }
        # Dictionary holding molar masses for reactants [g/mol]
        self.MMDict = {
            "LO2": 15.999,
            "N2O": 44.013,
            "ABS": 57.07,
            "Acrylic": 100.12
        }

        self.addCustomSpecies()
        self.C = CEA_Obj(oxName=oxName, fuelName=fuelName)



    def getOutput(self, Pc, OFRatio, ExpansionRatio, printOutput):
        output = self.C.get_full_cea_output(Pc=Pc, MR=OFRatio, eps=ExpansionRatio, short_output=0, output='siunits')
        if printOutput:
            print(output)

    def getChamberEquilibrium(self, Pc, OFRatio, ExpansionRatio, location):
        # gets all of the properties of the thingy (make sure not to ask for too much, you will upset him)
        # location: 0-injector (not supported) 1-chamber 2-throat 3-exit

        if location == 0:
            print("ERROR: properties at injector not supported")
        elif location == 1:
            self.Isp_Vac, self.C_star, self.T_flame, self.MW, self.gamma = self.C.get_IvacCstrTc_ChmMwGam(Pc=Pc, MR=OFRatio, eps=ExpansionRatio)  # get ISP [s], C* [ft/s], T [R], MW [g/mol], gamma [-]
            self.Cp, self.visc, self.thermCond, self.prandtl = self.C.get_Chamber_Transport(Pc, OFRatio, ExpansionRatio)  # get heat capacity [cal/g*K], viscosity[milliPoise], thermal conductivity [mcal/cm*s*K], Prandtl Number [-]

        elif location == 2:
            self.Isp_Vac, self.C_star, self.T_flame, self.MW, self.gamma = self.C.get_IvacCstrTc_ThtMwGam(Pc=Pc, MR=OFRatio, eps=ExpansionRatio)  # get ISP [s], C* [ft/s], T [R], MW [g/mol], gamma [-]
            self.Cp, self.visc, self.thermCond, self.prandtl = self.C.get_Throat_Transport(Pc, OFRatio, ExpansionRatio)  # get heat capacity [cal/g*K], viscosity[milliPoise], thermal conductivity [mcal/cm*s*K], Prandtl Number [-]
        elif location == 3:
            self.Isp_Vac, self.C_star, self.T_flame, self.MW, self.gamma = self.C.get_IvacCstrTc_exitMwGam(Pc=Pc, MR=OFRatio, eps=ExpansionRatio)  # get ISP [s], C* [ft/s], T [R], MW [g/mol], gamma [-]
            self.Cp, self.visc, self.thermCond, self.prandtl = self.C.get_Exit_Transport(Pc, OFRatio, ExpansionRatio)  # get heat capacity [cal/g*K], viscosity[milliPoise], thermal conductivity [mcal/cm*s*K], Prandtl Number [-]
        else:
            print("ERROR: location index out of bounds")
        self.T_flame = self.T_flame/1.8  # convert to Kelvin
        self.C_star = self.C_star / 3.28084  # convert characteristic velocity to [m/s]
        self.Cp = self.Cp * 4.186798  # get specific heat capacity [j/g*C]
        self.visc = self.visc/1000  # convert from milliPoise to Poise
        self.thermCond = self.thermCond * 418000  # convert to [W/m*K]

        return self.Isp_Vac, self.C_star, self.T_flame, self.MW, self.gamma, self.Cp, self.visc, self.thermCond, self.prandtl

    def getReactionHeat(self, Pc, OFRatio, ExpansionRatio, location):
        h_products = self.getProductsEnthalpy(Pc, OFRatio, ExpansionRatio, location)
        h_fuel = self.findHeatOfFormation(self.fuelName)/self.MMDict[self.fuelName]
        h_ox = self.findHeatOfFormation(self.oxName)/self.MMDict[self.oxName]

        X_fuel = 1/(OFRatio+1)  # fuel mass fraction
        X_ox = OFRatio/(OFRatio+1)  # oxidizer mass fraction

        h_reactants = X_fuel*h_fuel + X_ox*h_ox

        Q_total = h_products - h_reactants  # assuming adiabatic circumstances, Q_total = deltaH
        return Q_total

    def getProductsEnthalpy(self, Pc, OFRatio, ExpansionRatio, location):
        # calculates the enthalpy per gram of products
        # location: 0-injector 1-chamber 2-throat 3-exit

        molarMass, massFraction = self.C.get_SpeciesMassFractions(Pc=Pc, MR=OFRatio, eps=ExpansionRatio)

        # convert molarMass and massFraction into lists with species name in first column, and data in second
        molarMass = list(molarMass), list(molarMass.values())
        massFraction = list(massFraction), list(massFraction.values())

        # create an array containing the standard heat of formations of all species
        h_f = molarMass
        h_products = 0
        for i in range(len(molarMass[0])):
            h_f[1][i] = self.findHeatOfFormation(molarMass[0][i]) / molarMass[1][i]  # heat of formation [kJ/g]
            h_products += h_f[1][i] * massFraction[1][i][location]
        return h_products

    def findHeatOfFormation(self, species):
        try:
            return self.heatDict[species]
        except:
            print("ERROR:", species, "is not defined in heat of formation dictionary. Look in the init of the CEA class")

    def addCustomSpecies(self):
        # Add custom species to CEA
        card_str = """
                fuel ABS  C 3.85   H 4.85   N 0.43     wt%=100.00
                h,cal=14990    t(k)=298.15   rho=0.975
                """
        add_new_fuel('ABS', card_str)
        card_str = """
                fuel Acrylic  C 5   H 8   O 2     wt%=100.00
                h,cal=91396    t(k)=298.15   rho=1.18
                """
        add_new_fuel('Acrylic', card_str)