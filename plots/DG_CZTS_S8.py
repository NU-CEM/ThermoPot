################################################################################
#  Copyright Adam J. Jackson (2014)                                            #
#                                                                              #
#   This program is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   This program is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
################################################################################

def main():
    from materials import CZTS, Cu, Zn, Sn, S8
    import numpy as np
    T = np.linspace(100,1500,100)    # K
    P = np.array( np.logspace(1,7,100),ndmin=2).transpose() # Pa
    
    D_mu = CZTS.mu_kJ(T,P) - (2*Cu.mu_kJ(T,P) +
                                    Zn.mu_kJ(T,P) +
                                    Sn.mu_kJ(T,P) +
                                    0.5*S8.mu_kJ(T,P)
        )
    
    D_mu_label = '$\Delta G_f$ / kJ mol$^{-1}$'
    scale_range = [-380,-240]
    
    plot_potential(T,P,D_mu,D_mu_label,scale_range, filename='plots/DG_CZTS_S8.png')

def plot_potential(T,P,potential,potential_label,scale_range,filename=False,
                   precision="%d",T_units='K',P_units='Pa'):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rcParams['font.size'] = 16

    # Unit conversions (all calculations are in SI units, conversion needed for plots)
    if T_units=='K':
        x_values = T
        x_unitlabel = 'K'
    elif T_units=='C':
        x_values = T - 273.15
        x_unitlabel = '$^\circ$ C'
    else:
        raise ValueError('Invalid temperature unit: {0}'.format(T_units))


    if P_units=='Pa':
        y_values = P.flatten()
    elif P_units=='Bar' or P_units=='bar':
        y_values = P.flatten()*1E-5
    elif P_units=='mbar': 
        y_values = P.flatten()*1E-5*1E3
    elif P_units=='kPa':
        y_values = P.flatten()*1E-3
    elif P_units=='mmHg' or P_units=='torr':
        y_values = P.flatten()*760/(1.01325E5)
    else:
        raise ValueError('Invalid pressure unit: {0}.'.format(T_units))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    a = plt.contour(x_values,y_values,potential,10, linewidths = 1, colors = 'k')
    plt.pcolormesh(x_values,y_values,potential,
    #              cmap=plt.get_cmap('BuPu'),vmin=scale_range[0], vmax=scale_range[1])
                   cmap=plt.get_cmap('summer'),vmin=scale_range[0], vmax=scale_range[1])
    colours = plt.colorbar()
    plt.xlabel('Temperature / {0}'.format(x_unitlabel))    
    plt.ylabel('Pressure / {0}'.format(P_units))
    colours.set_label(potential_label, labelpad=20)
    ax.set_yscale('log')
    plt.clabel(a,fmt=precision)
    if filename:
        plt.savefig(filename,dpi=200)
    else:
        plt.show()

    
if __name__ == "__main__":
    main()
