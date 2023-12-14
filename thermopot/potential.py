import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
import numpy as np


class Potential:
    def __init__(self, potential, T, P):
        self.potential = potential
        self.T = T
        self.P = P

    def plot_TvsP(
        self,
        potential_label="$\Delta G_f$ (kJ mol$^{-1})$",
        scale_range=[-600, 0],
        filename=None,
        precision="%d",
        T_units="K",
        P_units="Pa",
        log_scale=True,
        sulphur_gas=False,
        gas_phase="S2",
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
        logscale determines if the y-axis Pressure is logarithmic
        """

        mpl.rcParams["font.family"] = "serif"
        mpl.rcParams["font.serif"] = "Times New Roman"
        mpl.rcParams["font.size"] = 20

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
        colormap = plt.get_cmap("GnBu")

        plt.pcolormesh(
            x_values,
            y_values,
            self.potential,
            cmap=colormap.reversed(),
            vmin=scale_range[0],
            vmax=scale_range[1],
            shading="auto",
            zorder=1,
        )
        colours = plt.colorbar()
        colours.set_label(potential_label, labelpad=20)

        if log_scale:
            ax.set_yscale("log")

        plt.xlabel("Temperature ({0})".format(x_unitlabel))
        plt.ylabel("Pressure ({0})".format(P_units))

        if sulphur_gas:
            T_tr_poly = [8.492e-01, 2.662e00, 3.849e01, 5.336e02]
            # print(self.T)
            pressure = self.P

            def T_tr(P):
                return np.polyval(T_tr_poly, np.log10(P))

            x = T_tr(pressure).flatten()

            if T_units == "C":
                x = x - 273.15

            if gas_phase == "S2":
                plt.fill_between(
                    x,
                    y_values,
                    9500000,
                    facecolor="w",
                    alpha=1,
                    zorder=4,
                    hatch="///",
                    linewidth=0,
                    edgecolor="0.8",
                )
                x1 = [
                    3,
                    419.15 - 265.15,
                ]  # LW: I have shifted this without fully understanding the logic of x1,y1,y2...
                y1 = [9500000, 9500000]
                y2 = [0.001, 0.001]
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
                plt.plot(x, y_values, "k--", linewidth=3)
                plt.xlim(min(x_values), max(x_values))
                plt.ylim(min(y_values), max(y_values))
                # plt.fill_between(x1, y1, y2, facecolor="none",edgecolor='k',hatch='/',zorder=5)
            elif gas_phase == "S8":
                plt.fill_between(
                    x,
                    y_values,
                    0.001,
                    facecolor="w",
                    alpha=1,
                    zorder=5,
                    hatch="///",
                    linewidth=0,
                    edgecolor="0.8",
                )
                x1 = [
                    950,
                    1000,
                ]  # LW: I have shifted this without fully understanding the logic of x1,y1,y2...
                y1 = [9500000, 9500000]
                y2 = [0.001, 0.001]
                plt.fill_between(
                    x1,
                    y1,
                    y2,
                    facecolor="w",
                    alpha=1,
                    zorder=3,
                    hatch="///",
                    linewidth=0,
                    edgecolor="0.8",
                )
                plt.plot(x, y_values, "k--", linewidth=3)
                plt.xlim(min(x_values), max(x_values))
                plt.ylim(min(y_values), max(y_values))
            elif gas_phase == "full":
                x1 = [
                    3,
                    1000,
                ]  # LW: I have shifted this without fully understanding the logic of x1,y1,y2...
                y1 = [9500000, 9500000]
                y2 = [0.001, 0.001]
                plt.fill_between(
                    x1,
                    y1,
                    y2,
                    facecolor="w",
                    alpha=1,
                    zorder=4,
                    linewidth=0,
                    edgecolor="0.8",
                )
                plt.plot(
                    x,
                    y_values,
                    color="k",
                    linestyle="--",
                    dashes=(7.6, 3.5),
                    linewidth=1.5,
                    zorder=5,
                )
                plt.xlim(min(x_values), max(x_values))
                plt.ylim(min(y_values), max(y_values))

        a = plt.contour(
            x_values,
            (y_values),
            self.potential,
            linewidths=0.7,
            colors="k",
            zorder=2,
        )
        plt.clabel(a, fmt=precision, zorder=3, fontsize=17)

        if filename:
            plt.savefig(filename, dpi=200, bbox_inches="tight")
        else:
            return plt
