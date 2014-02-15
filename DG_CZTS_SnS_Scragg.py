def main():
    from materials import CZTS, Cu2S_low, ZnS_zincblende, SnS, S2
    import numpy as np

    # Surface model from ab initio calcs and NIST data

    T = np.linspace(573.15,873.15,10)    # K
    P = np.array( np.logspace(-8,2,100),ndmin=2).transpose() # Pa
    
    D_mu = CZTS.mu_kJ(T,P) - (Cu2S_low.mu_kJ(T,P) +
                              ZnS_zincblende.mu_kJ(T,P) +
                              SnS.mu_kJ(T,P) +
                              0.5*S2.mu_kJ(T,P) )
    
    D_mu_label = '$\Delta G_f$ / kJ mol$^{-1}$'
    scale_range = [-50,70]

    # Stability lines from figure 5, J. J. Scragg et al., Chem. Mater. (2011) 23 4625-4633
    kinetic_data = np.genfromtxt('jscragg_2011.csv',delimiter=',',skip_header=1)
    # Convert log pressure to absolute pressure in mbar
    kinetic_data[:,1:] = np.power(10,kinetic_data[:,1:])
    # Columnwise conversion to SI units from deg C and mbar
    kinetic_data_si = (kinetic_data + [273.15, 0., 0.]) * [1., 100., 100.]
    
    plot_potential(T,P,D_mu,D_mu_label,scale_range, filename='plots/DG_CZTS_SnS_Scragg.png',
                   T_units='K', P_units='Pa', overlay=kinetic_data_si)

def plot_potential(T,P,potential,potential_label,scale_range,filename=False,
                   precision="%d",T_units='K',P_units='Pa',overlay=False):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rcParams['font.size'] = 16

    # Unit conversions (all calculations are in SI units, conversion needed for plots)
    x_values = si_to_other(T,T_units)
    y_values = si_to_other(P.flatten(),P_units)
    if T_units=='K':
        x_unitlabel = 'K'
    elif T_units=='C':
        x_unitlabel = '$^\circ$ C'
    else:
        raise ValueError('Invalid temperature unit: {0}'.format(T_units))

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
    
    # Overlay additional data lines
    if overlay.any():
        #        ax.plot(overlay[:,0],overlay[:,1:],'r')
        ax.fill_between(si_to_other(overlay[:,0],T_units),si_to_other(overlay[:,1],P_units),
                        si_to_other(overlay[:,2],P_units),edgecolor='#444444',facecolor='none',hatch='xx')
        ax.axis([min(x_values),max(x_values),min(y_values),max(y_values)])
        

    if filename:
        plt.savefig(filename,dpi=200)
    else:
        plt.show()

def si_to_other(data, units):
    # Pressure conversions
    if units=='Pa':
        pass
    elif units=='Bar' or units=='bar':
        data = data*1E-5
    elif units=='mbar': 
        data = data*1E-5*1E3
    elif units=='kPa':
        data = data*1E-3
    elif units=='mmHg' or units=='torr':
        data = data*760/(1.01325E5)
    # Temperature conversions
    elif units=='K':
        pass
    elif units=='C':
        data = data - 273.15
    else:
        raise ValueError('Invalid unit: {0}.'.format(units))

    return data    

if __name__ == "__main__":
    main()
