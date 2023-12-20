import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import numpy as np

class Potentials:
    def __init__(self, *potentials):
        self.potentials = potentials
        self.minimum_potential = self.find_potential_minimum()
        self.T = self.potentials[0].T
        self.P = self.potentials[0].P

    # TODO: check that all the potentials are plottable over T/P range

    def plot_TvsP(
        self,
        material_labels=None,
        T_units="K",
        P_units="Pa",
        log_scale=True,
    ):
        """
        T is an array e.g. np.linspace(100, 1500, 100)  # K
        P is an array orthogonal to T. e.g. np.array(np.logspace(1, 7, 100), ndmin=2).transpose()  # Pa
        potential is returned from a reactions.reaction method called for an instance with attributes T,P.
        If T has length m and P has length n, P will be a 2D array with dimensions m x n.
        e.g. reactions.reaction({Ba:1,S:2}, {BaS2:1}},temperature=T,pressure=P).Dmu_eV_pbesol()
        potential_label is the label of the contour colorbar e.g. '$\Delta G_f$ / kJ mol$^{-1}$'
        scale_range is the scale of the colorbar e.g. [-380, -240]
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
            x_unitlabel = "$^\circ$ C"
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
        )
        # TODO: sort the colour map out so consistent with grid. Now ranges from 0 to 1

        # Set borders in the interval [0, 1]
        bound = np.linspace(0, 1, len(material_labels))

        # AT THE MOMENT THIS IS BROKEN!!!!
        # plt.legend(
        #    [mpatches.Patch(color=colormap(i)) for i in bound],
        #    ["{:s}".format(material_labels[i]) for i in range(len(material_labels))],
        # )

        plt.xlabel("Temperature / {0}".format(x_unitlabel))
        plt.ylabel("Pressure / {0}".format(P_units))
        if log_scale:
            ax.set_yscale("log")

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
