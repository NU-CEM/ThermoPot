from thermopot import potential


class Reaction:
    """
    Class for reaction data

    Sets properties:
    -------------------
    reaction.reactants     (Dict relating reactant materials to a number of formula units)
    reaction.products      (Dict relating product materials to a number of formular units)

    Sets methods:
    -------------------
    reaction.DH(T,P) : Enthalpy of formation
    reaction.DU(T,P) : Internal energy change
    reaction.Dmu(T,P) : Gibbs free energy of formation
    """

    def __init__(
        self,
        reactants_dictionary,
        products_dictionary,
        temperature=298.15,
        pressure=1e5,
        fu=1,
    ):
        """
        reactants_dictionary and products dictionary takes the form { class_instance : formula units }
        and can have an arbitrary number of key-value pairs. `Class instance` is an instance of the `materials.solid`
        or `materials.ideal_gas` classes.

        temperature is provided in kelvin, pressure is provided in Pa.

        fu is the number of formula units of the final reactant(s). It is used to scale the calculated changes in energy/enthalpy.

        """
        self.reactants = reactants_dictionary
        self.products = products_dictionary
        self.T = temperature
        self.P = pressure
        self.fu = fu

    def DH(self, T=None, P=None, xc="pbesol", units="eV"):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_enthalpy, products_enthalpy = 0, 0
        for material, fu in self.reactants.items():
            reactants_enthalpy += material.H(T, P, xc=xc, units=units) * fu
        for material, fu in self.products.items():
            products_enthalpy += material.H(T, P, xc=xc, units=units) * fu

        return potential.Potential(
            (products_enthalpy - reactants_enthalpy) / self.fu, T, P
        )

    def DU(self, T=None, P=None, xc="pbesol", units="eV"):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.U(T, xc=xc, units=units) * fu
        for material, fu in self.products.items():
            products_energy += material.U(T, xc=xc, units=units) * fu

        return potential.Potential((products_energy - reactants_energy) / self.fu, T, P)

    def Dmu(self, T=None, P=None, xc="pbesol", units="eV"):

        T = self.T if T is None else T
        P = self.P if P is None else P

        reactants_energy, products_energy = 0, 0
        for material, fu in self.reactants.items():
            reactants_energy += material.mu(T, P, xc=xc, units=units) * fu
        for material, fu in self.products.items():
            products_energy += material.mu(T, P, xc=xc, units=units) * fu

        return potential.Potential((products_energy - reactants_energy) / self.fu, T, P)
