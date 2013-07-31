"""Calculate enthalpies of formation for CZTS over a range of temperatures at standard pressure."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import constants # Library provides 2010 CODATA values for physical constants
from interpolate_thermal_property import get_potential_aims, get_potential_nist_table
from aims import pbesol_energy_eV, fu_cell  # Data file of FHI-aims computed systems and results
from aims import volume as cell_volume

def main():
########## Alias physical constants ###########
    eV2J = constants.physical_constants['electron volt-joule relationship'][0]
    k = constants.physical_constants['Boltzmann constant in eV/K'][0]

################ Set conditions ###############
    T = np.linspace(100,1000,num=100) # K
    P = np.array(10**np.linspace(4,6,num=50),ndmin=2).transpose()    # Pa
############# Data and parameters from DFT calcs ################
 
    ### Copy energies and convert to J: ###
    pbesol_energy_J = pbesol_energy_eV.copy()
    for species in pbesol_energy_eV:
        pbesol_energy_J[species] = pbesol_energy_eV[species]*eV2J*constants.N_A
    ### thermodynamic data for periodic solids ###
    thermo_file=dict()
    U_func = dict()
    H_func = dict()
    TS_func = dict()
    S_func = dict()

    phonon_species = ('czts','Cu','Zn','Sn','alpha_S')
    for species in phonon_species:
        thermo_file[species] = 'phonopy_output/' + species + '.dat'
        U_func[species] = get_potential_aims(thermo_file[species],'U')
        TS_func[species] = get_potential_aims(thermo_file[species],'TS')

    ### thermodynamic data for non-periodic components ###
    gas_species = ['S8']
    for species in gas_species:
        thermo_file[species] = 'nist_janaf/' + species + '.dat'
        H_func[species] = get_potential_nist_table(thermo_file[species],'H')
        S_func[species] = get_potential_nist_table(thermo_file[species],'S')
        
############## Enthalpy for each species #####################
    H = dict()
    mu = dict()
    for species in phonon_species:
        H[species] = pbesol_energy_eV[species] + U_func[species](T) + (P*cell_volume[species]*1E-30)/(eV2J*constants.N_A)
        mu[species] = H[species] - TS_func[species](T)
    for species in gas_species:
        H[species] = pbesol_energy_eV[species] + (H_func[species](T))/(eV2J*constants.N_A) # Assuming ideality
        mu[species] = H[species] - T*S_func[species](T)/(eV2J*constants.N_A) + k*T*np.log(P/1E5)
######### Calculate some formation energies! #################
    formation_S8_eV = (pbesol_energy_eV['czts']/fu_cell['czts'] -
                          (2*pbesol_energy_eV['Cu']/fu_cell['Cu'] +
                           pbesol_energy_eV['Zn']/fu_cell['Zn'] +
                           pbesol_energy_eV['Sn']/fu_cell['Sn'] +
                           0.5*pbesol_energy_eV['S8']
                           ))
    print formation_S8_eV
    DH_S8_eV = (H['czts']/fu_cell['czts'] -
                (2*H['Cu']/fu_cell['Cu'] +
                 H['Zn']/fu_cell['Zn'] +
                 H['Sn']/fu_cell['Sn'] +
                 0.5*H['S8']
                 ))
    print DH_S8_eV
    DG_S8_eV = (mu['czts']/fu_cell['czts'] -
                (2*mu['Cu']/fu_cell['Cu'] +
                 mu['Zn']/fu_cell['Zn'] +
                 mu['Sn']/fu_cell['Sn'] +
                 0.5*mu['S8']
                 ))
    print DG_S8_eV
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    a = plt.contour(T,P.flatten(),DG_S8_eV,10, linewidths = 0.5, colors = 'k')
    plt.pcolormesh(T,P.flatten(),DG_S8_eV,cmap=plt.get_cmap('BuPu'))
    colours = plt.colorbar()
    plt.xlabel('Temperature / K')
    plt.ylabel('Pressure / Pa')
    plt.title('2 Cu + Zn + Sn + 0.5 S8 --> CZTS')
    colours.set_label('DG / eV')
    ax.set_yscale('log')
    plt.clabel(a)
    plt.show()
    plt.savefig("formation_from_S8.eps")
    return 0

if __name__=="__main__":
    main()



