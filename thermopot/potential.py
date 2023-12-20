import matplotlib.pyplot as plt
import matplotlib as mpl


class Potential:
    def __init__(self, potential, T, P):
        self.potential = potential
        self.T = T
        self.P = P

    def plot_TvsP(
        self,
        potential_label="$\Delta G_f$ / kJ mol$^{-1}$",
        scale_range=[-600, 0],
        precision="%d",
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
        logscale determines if the y-axis Pressure is logarithmic
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

        a = plt.contour(
            x_values, y_values, self.potential, 10, linewidths=1, colors="k"
        )
        plt.pcolormesh(
            x_values,
            y_values,
            self.potential,
            cmap=colormap,
            vmin=scale_range[0],
            vmax=scale_range[1],
            shading="auto",
        )
        colours = plt.colorbar()
        colours.set_label(potential_label, labelpad=20)

        if log_scale:
            ax.set_yscale("log")
        plt.clabel(a, fmt=precision)

        plt.xlabel("Temperature / {0}".format(x_unitlabel))
        plt.ylabel("Pressure / {0}".format(P_units))

        return plt
