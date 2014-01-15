def main():
    from materials import CZTS, Cu, Zn, Sn, S8
    import numpy as np
    T = np.linspace(100,1000,100)    # K
    P = np.array( np.logspace(1,7,100),ndmin=2).transpose() # Pa
    
    D_mu = CZTS.mu_kJ(T,P) - (2*Cu.mu_kJ(T,P) +
                                    Zn.mu_kJ(T,P) +
                                    Sn.mu_kJ(T,P) +
                                    0.5*S8.mu_kJ(T,P)
        )
    
    D_mu_label = '$\Delta G_f$ / kJ mol$^{-1}$'
    scale_range = [-380,-240]
    
    plot_potential(T,P,D_mu,D_mu_label,scale_range)

def plot_potential(T,P,potential,potential_label,scale_range):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = ['Charter']
    mpl.rcParams['font.size'] = 16
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    a = plt.contour(T,P.flatten(),potential,10, linewidths = 1, colors = 'k')
    plt.pcolormesh(T,P.flatten(),potential,
                   cmap=plt.get_cmap('BuPu'),vmin=scale_range[0], vmax=scale_range[1])
    colours = plt.colorbar()
    plt.xlabel('Temperature / K')
    plt.ylabel('Pressure / Pa')
    colours.set_label(potential_label)
    ax.set_yscale('log')
    plt.clabel(a,fmt="%d<")
    
    plt.show()
    
if __name__ == "__main__":
    main()
