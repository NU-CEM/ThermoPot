import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import numpy as np
import scipy


class Potentials:
    def __init__(self, *potentials):
        self.potentials = potentials
        self.minimum_potential = self.find_potential_minimum()
        self.T = self.potentials[0].T
        self.P = self.potentials[0].P

    # TODO: check that all the potentials are plottable over T/P range
    # TODO: check that all the potentials are calculated for the same T/P range.

    def plot_TvsP(
        self,
        material_labels=None,
        filename=False,
        T_units="K",
        P_units="Pa",
        log_scale=True,
        sulphur_gas=False,
        gas_phase="S2",
        melting_point=False,
    ):
        """
        T is an array e.g. np.linspace(100, 1500, 100)  # K
        P is an array orthogonal to T. e.g. np.array(np.logspace(1, 7, 100), ndmin=2).transpose()  # Pa
        potential is returned from a reactions.reaction method called for an instance with attributes T,P.
        If T has length m and P has length n, P will be a 2D array with dimensions m x n.
        e.g. reactions.reaction({Ba:1,S:2}, {BaS2:1}},temperature=T,pressure=P).Dmu_eV_pbesol()
        potential_label is the label of the contour colorbar e.g. '$\Delta G_f$ / kJ mol$^{-1}$'
        scale_range is the scale of the colorbar e.g. [-380, -240]
        filename is the output filename e.g. 'plots/Dmu-BaS2-Ba-S2.png'. If not provided `plt.show()` is called.
        log_scale determines if the Pressure y-axis is logarithmic
        """

        mpl.rcParams["font.family"] = "serif"
        mpl.rcParams["font.serif"] = "Times New Roman"
        mpl.rcParams["font.size"] = 16

        # Unit conversions (all calculations are in SI units, conversion needed for plots)

        if T_units == "K":
            x_values = self.T
            x_unitlabel = "K"
        elif T_units == "C":
            x_values = self.T - 273.15
            x_unitlabel = "$\degree$C"
        else:
            raise ValueError("Invalid temperature unit: {0}".format(T_units))

        if P_units == "Pa":
            y_values = self.P.flatten()
        elif P_units == "Bar" or P_units == "bar":
            y_values = self.P.flatten() * 1e-5
        elif P_units == "mbar":
            y_values = self.P.flatten() * 1e-5 * 1e3
        elif P_units == "kPa":
            y_values = self.P.flatten() * 1e-3
        elif P_units == "mmHg" or P_units == "torr":
            y_values = self.P.flatten() * 760 / (1.01325e5)
        else:
            raise ValueError("Invalid pressure unit: {0}.".format(T_units))

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        colormap = plt.get_cmap("summer")

        potential = self.find_potential_minimum()
        plt.pcolormesh(
            x_values,
            y_values,
            potential / (len(material_labels) - 1),
            cmap=colormap,
            shading="auto",
            zorder=1,
        )

        # AT THE MOMENT THIS IS BROKEN!!!!
        # plt.legend([mpatches.Patch(color=colormap(i)) for i in bound],
        #    ["{:s}".format(material_labels[i]) for i in range(len(material_labels))],
        # )

        plt.xlabel("Temperature ({0})".format(x_unitlabel))
        plt.ylabel("Pressure ({0})".format(P_units))

        if log_scale:
            ax.set_yscale("log")

        if filename:
            plt.savefig(filename, dpi=200)
        else:
            return plt

    def find_potential_minimum(self):
        assert (
            len(set([potential.potential.shape for potential in self.potentials])) == 1
        ), "potential arrays must have the same dimension"

        minimum_potential = self.potentials[0].potential
        for i, potential in enumerate(self.potentials):
            minimum_potential = np.minimum(
                minimum_potential, self.potentials[i + 1].potential
            )
            if i + 2 == len(self.potentials):
                break

        for i, potential in enumerate(self.potentials):
            minimum_potential[potential.potential == minimum_potential] = i

        return minimum_potential

    def polynomial_intersection(
        self,
        ratio=4,
        pressures=np.logspace(0,7,100),
        order=4,
        ref_pressure_for_fit=1E5
        ):
        """
        Conducts polynomial fittings to two potential datasets.
        Calculates the temperature(s) at which the two polynomial fits intersects for given pressure(s).
        Fits another polynomial to the transition-temp vs pressure data. 
        This can then be evaluated and overlaid on potential contour plots (not implemented here).
        
        pressures is given a iterable with units in pascals. It defaults to np.logspace(0,7,100).
        ratio is num_atoms_potential_1 / num_atoms_potential_2. It defaults to 4.
        order is the order of the polynomial fitting. It defaults to 4.
        ref_pressure_for_fit is a float and is given in pascals. It defaults to one bar.

        If the calculated potential is not available at this pressure then an error is raised.

        This broadly follows code published by Adam Jackson: 
        https://github.com/WMD-group/sulfur-model/blob/master/parameterisation.ipynb
        """

        try:
            len(self.potentials) != 2
        except: 
            print("this method can only be called for a Potentials object storing two potentials")

        # generate polynomial fits to calculated data
        poly1 = self.potentials[0].polyfit_potential(order,ref_pressure_for_fit)
        poly2 = self.potentials[1].polyfit_potential(order,ref_pressure_for_fit)
        
        # xp1 is coefficient for potential 1
        # xp2 is coefficient for potential 2
        
        def build_poly_to_solve(poly1, poly2, ratio, pressure):
            
            k = scipy.constants.physical_constants['Boltzmann constant in eV/K'][0]
            
            return ([xp1-ratio*xp2 for xp1,xp2 in zip(poly1[:-2],poly2[:-2])]
                           + [(poly1[-2]-ratio*poly2[-2]-(ratio-1)*k*np.log(pressure/ref_pressure_for_fit))]
                           + [poly1[-1] - ratio*poly2[-1]])

        # list to store temperatures at which intersect
        T_intersect = []

        # if pressure is float then convert to iterable
        if not hasattr(pressures, '__iter__'):
            pressures = [pressures]
        
        for pressure in pressures:
            poly_to_solve = build_poly_to_solve(poly1, poly2, ratio, pressure)
            roots = np.roots(poly_to_solve)
            
            best_root = False
            for root in roots:
                if root == root.real and root > 0 and best_root == False:
                    best_root = root.real
                elif root == root.real and root > 0 and root < best_root:
                    best_root = root.real
            if best_root:
                T_intersect.append(best_root)
            else:
                T_intersect.append('NaN')
                
        T_intersect = [t.real for t in T_intersect]
        
        T_intersect_fit = np.polyfit(np.log10(pressures), T_intersect, 3)
        
        return [T_intersect_fit,T_intersect]

        