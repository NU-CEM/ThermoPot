import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_TvsP(T, P, *potential, potential_label="$\Delta G_f$ / kJ mol$^{-1}$", scale_range=[-600, 0], filename=False,
                   precision="%d", T_units='K', P_units='Pa'):
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

    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rcParams['font.size'] = 16

    # Unit conversions (all calculations are in SI units, conversion needed for plots)
    if T_units == 'K':
        x_values = T
        x_unitlabel = 'K'
    elif T_units == 'C':
        x_values = T - 273.15
        x_unitlabel = '$^\circ$ C'
    else:
        raise ValueError('Invalid temperature unit: {0}'.format(T_units))

    if P_units == 'Pa':
        y_values = P.flatten()
    elif P_units == 'Bar' or P_units == 'bar':
        y_values = P.flatten() * 1E-5
    elif P_units == 'mbar':
        y_values = P.flatten() * 1E-5 * 1E3
    elif P_units == 'kPa':
        y_values = P.flatten() * 1E-3
    elif P_units == 'mmHg' or P_units == 'torr':
        y_values = P.flatten() * 760 / (1.01325E5)
    else:
        raise ValueError('Invalid pressure unit: {0}.'.format(T_units))

    if len(potential) == 1:
        potential = potential

    else:
        find_potentials_minimum(potential)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    a = plt.contour(x_values, y_values, potential, 10, linewidths=1, colors='k')
    plt.pcolormesh(x_values, y_values, potential,
                   cmap=plt.get_cmap('summer'), vmin=scale_range[0], vmax=scale_range[1])
    colours = plt.colorbar()
    plt.xlabel('Temperature / {0}'.format(x_unitlabel))
    plt.ylabel('Pressure / {0}'.format(P_units))
    colours.set_label(potential_label, labelpad=20)
    ax.set_yscale('log')
    plt.clabel(a, fmt=precision)
    if filename:
        plt.savefig(filename, dpi=200)
    else:
        plt.show()

def find_potentials_minimum(*potential):

    assert (len(set([array.shape() for array in potentials])) == 1), "potential arrays must have the same dimension"
    pass
    # need to think about this more. Need to label potentials for a key object.