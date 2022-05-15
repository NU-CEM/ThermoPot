import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

class Potentials():

    def __init__(self,*potentials):

        self.potentials = potentials
        self.minimum_potential = self.find_potential_minimum()

    # TODO: check that all the potentials are plottable over T/P range

    def plot_TvsP(self,
                T = None,
                P = None,
                material_labels=None,
                filename=False,
                T_units="K",
                P_units="Pa"
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
        """

        mpl.rcParams["font.family"] = "serif"
        mpl.rcParams["font.serif"] = "Times New Roman"
        mpl.rcParams["font.size"] = 16

        # Unit conversions (all calculations are in SI units, conversion needed for plots)
        if T_units == "K":
            x_values = T
            x_unitlabel = "K"
        elif T_units == "C":
            x_values = T - 273.15
            x_unitlabel = "$^\circ$ C"
        else:
            raise ValueError(
                "Invalid temperature unit: {0}".format(T_units))

        if P_units == "Pa":
            y_values = P.flatten()
        elif P_units == "Bar" or P_units == "bar":
            y_values = P.flatten() * 1e-5
        elif P_units == "mbar":
            y_values = P.flatten() * 1e-5 * 1e3
        elif P_units == "kPa":
            y_values = P.flatten() * 1e-3
        elif P_units == "mmHg" or P_units == "torr":
            y_values = P.flatten() * 760 / (1.01325e5)
        else:
            raise ValueError("Invalid pressure unit: {0}.".format(T_units))

        fig = plt.figure()
        colormap = plt.get_cmap("summer")

        potential = self.find_potential_minimum()
        plt.pcolormesh(
            x_values,
            y_values,
            potential,
            cmap=colormap
        )
        plt.legend(
            [mpl.patches.Patch(color=colormap(b)) for b in range(len(
                potential))],
            material_labels)

        plt.xlabel("Temperature / {0}".format(x_unitlabel))
        plt.ylabel("Pressure / {0}".format(P_units))

        if filename:
            plt.savefig(filename, dpi=200)
        else:
            plt.show()


    def find_potential_minimum(self):

        assert (
                len(set([array.shape for array in self.potentials])) == 1
        ), "potential arrays must have the same dimension"

        minimum_potential = self.potentials[0]
        for i, potential in enumerate(self.potentials):
            minimum_potential = np.minimum(minimum_potential,
                                           self.potentials[i + 1])
            if i + 2 == len(self.potentials):
                break

        for i, potential in enumerate(self.potentials):
            minimum_potential[potential == minimum_potential] = i

        return minimum_potential


