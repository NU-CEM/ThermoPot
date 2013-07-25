"""Calculate enthalpies of formation for CZTS over a range of temperatures at standard pressure."""
import numpy as np
from scipy import constants # Library provides 2010 CODATA values for physical constants
from interpolate_thermal_property import get_potential_aims, get_potential_nist_table
from aims import pbesol_energy_eV, fu_cell  # Data file of FHI-aims computed systems and results

def main():
########## Alias physical constants ###########
    eV2J = constants.physical_constants['electron volt-joule relationship'][0]
    k = constants.physical_constants['Boltzmann constant in eV/K'][0]

################ Set conditions ###############
    T = 298.15 # K
    P = 1e5    # Pa
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
    V = 1E-30 # Need to import correct cell volume from geometry file
    H = dict()
    mu = dict()
    for species in phonon_species:
        H[species] = pbesol_energy_eV[species] + U_func[species](T) + (P*V)/(eV2J*constants.N_A)
        mu[species] = H[species] - TS_func[species](T)
    for species in gas_species:
        H[species] = pbesol_energy_eV[species] + (H_func[species](T))/(eV2J*constants.N_A) # Assuming ideality
        mu[species] = H[species] - T*S_func[species](T) + k*T*np.log(P/1E5)
######### Calculate some formation energies! #################
    formation_alpha_eV = (pbesol_energy_eV['czts']/fu_cell['czts'] -
                          (2*pbesol_energy_eV['Cu']/fu_cell['Cu'] +
                           pbesol_energy_eV['Zn']/fu_cell['Zn'] +
                           pbesol_energy_eV['Sn']/fu_cell['Sn'] +
                           4*pbesol_energy_eV['alpha_S']/fu_cell['alpha_S']
                           ))
    
    DH_alpha_eV = (H['czts']/fu_cell['czts'] -
                (2*H['Cu']/fu_cell['Cu'] +
                 H['Zn']/fu_cell['Zn'] +
                 H['Sn']/fu_cell['Sn'] +
                 4*H['alpha_S']/fu_cell['alpha_S']
                 ))
    DG_alpha_eV = (mu['czts']/fu_cell['czts'] -
                (2*mu['Cu']/fu_cell['Cu'] +
                 mu['Zn']/fu_cell['Zn'] +
                 mu['Sn']/fu_cell['Sn'] +
                 4*mu['alpha_S']/fu_cell['alpha_S']
                 ))


    ###### Formatted print: heading with underline followed by values in eV and kJ #####

    alphastring = "2 Cu + Zn + Sn + 4 S(alpha) -> CZTS"
    print alphastring + "\n" + len(alphastring) * "="
    print "Crude DFT energy @ {0:3.5g} K, {1:3.3g} bar:".format(T,P/1E5).ljust(50)+"{0:5.3g} eV".format(formation_alpha_eV)
    print "{0}{1:5.3g} kJ/mol".format((50*" "),formation_alpha_eV*constants.N_A*eV2J/1000)
    print "Enthalpy change @ {0:3.5g} K, {1:3.3g} bar:".format(T,P/1E5).ljust(50)+"{0:5.3g} eV".format(DH_alpha_eV)
    print "{0}{1:5.3g} kJ/mol".format((50*" "),DH_alpha_eV*constants.N_A*eV2J/1000)
    print "Gibbs free energy @ {0:3.5g} K, {1:3.3g} bar:".format(T,P/1E5).ljust(50)+"{0:5.3g} eV".format(DG_alpha_eV)
    print "{0}{1:5.3g} kJ/mol".format((50*" "),DG_alpha_eV*constants.N_A*eV2J/1000)
                    
    return 0




if __name__=="__main__":
    main()
