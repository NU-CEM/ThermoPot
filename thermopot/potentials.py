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

        if sulphur_gas:
            T_tr_poly = [8.492e-01, 2.662e00, 3.849e01, 5.336e02]
            # print(self.T)
            pressure = self.P

            def T_tr(P):
                return np.polyval(T_tr_poly, np.log10(P))

            x = T_tr(pressure).flatten()
            
            if T_units == "C":
                x = x - 273.15

            plt.plot(x, y_values, "k--", linewidth=3)
            plt.xlim(min(x_values), max(x_values))
            if gas_phase == "S2":
                plt.fill_between(
                    x,
                    y_values,
                    10000000,
                    facecolor="w",
                    alpha=1,
                    zorder=4,
                    hatch="///",
                    linewidth=0,
                    edgecolor="0.8",
                )
                x1 = [1, 419.1596 - 273.15]
                y1 = [9999900, 9999900]
                y2 = [0.0001, 0.0001]
                plt.fill_between(
                    x1,
                    y1,
                    y2,
                    facecolor="w",
                    alpha=1,
                    zorder=4,
                    hatch="///",
                    linewidth=0,
                    edgecolor="0.8",
                )
                # plt.fill_between(x1, y1, y2, facecolor="none",edgecolor='k',hatch='/',zorder=5)
            if gas_phase == "S8":
                plt.fill_between(
                    x,
                    y_values,
                    0.001,
                    facecolor="w",
                    alpha=1,
                    zorder=4,
                    hatch="///",
                    linewidth=0,
                    edgecolor="0.8",
                )

            # resolution = 1000
            # temp = np.linspace(self.T[0], self.T[-1], resolution)
            # plt.plot(T_tr(y_values),(y_values),'k--', linewidth=3)

        #        if sulphur_gas:
        #            T_tr_poly = [8.492e-01, 2.662e+00, 3.849e+01, 5.336e+02]
        #
        #            def T_tr(P):
        #                return np.polyval(T_tr_poly, np.log10(self.P.flatten()))
        #
        #            resolution = 1000
        #            temp = np.linspace(self.T[0], self.T[-1], resolution)
        #            plt.plot(T_tr(y_values),(y_values),'k--', linewidth=3)
        #        # TODO: sort the colour map out so consistent with grid. Now ranges from 0 to 1

        # Set borders in the interval [0, 1]
        # bound = np.linspace(0, 1, len(material_labels))

        # AT THE MOMENT THIS IS BROKEN!!!!
        # plt.legend([mpatches.Patch(color=colormap(i)) for i in bound],
        #    ["{:s}".format(material_labels[i]) for i in range(len(material_labels))],
        # )

        plt.xlabel("Temperature ({0})".format(x_unitlabel))
        plt.ylabel("Pressure ({0})".format(P_units))
        plt.ylim([0.001, 10000000])
        plt.xlim([0, 1000])
        if melting_point == True:
            plt.axvline(x=554, zorder=20, color="k")

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
