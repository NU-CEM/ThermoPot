def main():
    from materials import CZTS_kesterite, CZTS_stannite
    import numpy as np

    from DG_CZTS_S8 import plot_potential

    T = np.linspace(100,1500,100)    # K
    P = np.array( np.logspace(1,7,100),ndmin=2).transpose() # Pa
    
    D_mu = CZTS_stannite.mu_kJ(T,P) - CZTS_kesterite.mu_kJ(T,P)
    
    D_mu_label = '$\Delta G_f$ / kJ mol$^{-1}$'
    scale_range = [2,4]
    
    plot_potential(T,P,D_mu,D_mu_label,scale_range, filename='plots/DG_stannite.png', precision="%.1f")

if __name__ == "__main__":
    main()
