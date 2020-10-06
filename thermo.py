# Thermodynamics and chemical equilibrium
# Written by: Thomas Kulin
# October 2, 2020

from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant


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