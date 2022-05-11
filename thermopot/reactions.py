class Reaction():
    """
    Class for reaction data

    Sets properties:
    -------------------
    reaction.reactants     (Dict relating reactant materials to a number of formula units)
    reaction.products      (Dict relating product materials to a number of formular units)

    Sets methods:
    -------------------
    reaction.DH_eV_pbesol(T,P), reaction.DH_kJ_pbesol(T,P), reaction.DH_eV_hse06(T,P), reaction.DH_kJ_hse06(T,P) : Enthalpy of formation
    reaction.DU_eV_pbesol(T,P), reaction.DU_kJ_pbesol(T,P), reaction.DU_eV_hse06(T,P), reaction.DU_kJ_hse06(T,P) : Internal energy change
    reaction.Dmu_eV_pbesol(T,P), reaction.Dmu_kJ_pbesol(T,P), reaction.Dmu_eV_hse06(T,P), reaction.Dmu_kJ_hse06(T,P) : Gibbs free energy of formation
    """
    def __init__(self, reactants_dictionary, products_dictionary,temperature=298.15,pressure=1E5):
        """
        reactants_dictionary and products dictionary takes the form { class_instance : formula units }
        and can have an arbitrary number of key-value pairs. `Class instance` is an instance of the `materials.solid`
        or `materials.ideal_gas` classes.

        temperature is provided in kelvin, pressure is provided in Pa.

        """
        self.reactants = reactants_dictionary
        self.products = products_dictionary
        self.T = temperature
        self.P = pressure

    def DH_eV_pbesol(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_enthalpy, products_enthalpy = 0, 0
        for material, fu in self.reactants.items():
            reactants_enthalpy += material.H_eV(T,P,xc="pbesol")*fu
        for material, fu in self.products.items():
            products_enthalpy += material.H_eV(T,P,xc="pbesol")*fu

        return products_enthalpy - reactants_enthalpy

    def DH_kJ_pbesol(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_enthalpy, products_enthalpy = 0, 0
        for material, fu in self.reactants.items():
            reactants_enthalpy += material.H_kJ(T, P, xc="pbesol") * fu
        for material, fu in self.products.items():
            products_enthalpy += material.H_kJ(T, P, xc="pbesol") * fu

        return products_enthalpy - reactants_enthalpy

    def DH_eV_hse06(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_enthalpy, products_enthalpy = 0, 0
        for material, fu in self.reactants.items():
            reactants_enthalpy += material.H_eV(T,P,xc="hse06")*fu
        for material, fu in self.products.items():
            products_enthalpy += material.H_eV(T,P,xc="hse06")*fu

        return products_enthalpy - reactants_enthalpy

    def DH_kJ_hse06(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_enthalpy, products_enthalpy = 0, 0
        for material, fu in self.reactants.items():
            reactants_enthalpy += material.H_kJ(T, P, xc="hse06") * fu
        for material, fu in self.products.items():
            products_enthalpy += material.H_kJ(T, P, xc="hse06") * fu

        return products_enthalpy - reactants_enthalpy

    def DU_eV_pbesol(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.U_eV(T, P, xc="pbesol") * fu
        for material, fu in self.products.items():
            products_energy += material.U_eV(T, P, xc="pbesol") * fu

        return products_energy - reactants_energy

    def DU_kJ_pbesol(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.U_kJ(T, P, xc="pbesol") * fu
        for material, fu in self.products.items():
            products_energy += material.U_kJ(T, P, xc="pbesol") * fu

        return products_energy - reactants_energy

    def DU_eV_hse06(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.U_eV(T, P, xc="hse06") * fu
        for material, fu in self.products.items():
            products_energy += material.U_eV(T, P, xc="hse06") * fu

        return products_energy - reactants_energy

    def DU_kJ_hse06(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.U_kJ(T, P, xc="hse06") * fu
        for material, fu in self.products.items():
            products_energy += material.U_kJ(T, P, xc="hse06") * fu

        return products_energy - reactants_energy

    def Dmu_eV_pbesol(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.mu_eV(T, P, xc="pbesol") * fu
        for material, fu in self.products.items():
            products_energy += material.mu_eV(T, P, xc="pbesol") * fu

        return products_energy - reactants_energy

    def Dmu_kJ_pbesol(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.mu_kJ(T, P, xc="pbesol") * fu
        for material, fu in self.products.items():
            products_energy += material.mu_kJ(T, P, xc="pbesol") * fu

        return products_energy - reactants_energy

    def Dmu_eV_hse06(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.mu_eV(T, P, xc="hse06") * fu
        for material, fu in self.products.items():
            products_energy += material.mu_eV(T, P, xc="hse06") * fu

        return products_energy - reactants_energy

    def Dmu_kJ_hse06(self, T=None, P=None):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.mu_kJ(T, P, xc="hse06") * fu
        for material, fu in self.products.items():
            products_energy += material.mu_kJ(T, P, xc="hse06") * fu

        return products_energy - reactants_energy


