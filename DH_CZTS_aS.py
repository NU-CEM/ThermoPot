"""Calculate enthalpies of formation for CZTS over a range of temperatures at standard pressure."""
from scipy import constants # Library provides 2010 CODATA values for physical constants
from interpolate_thermal_property import get_potential_aims
from aims import pbesol_energy_eV, fu_cell  # Data file of FHI-aims computed systems and results

def main():
################ Set conditions ###############
    T = 298.15 # K
    P = 1e5    # Pa
############# Data and parameters from DFT calcs ################
 
    ### Copy energies and convert to J: ###
    eV2J = constants.physical_constants['electron volt-joule relationship'][0]
    pbesol_energy_J = pbesol_energy_eV.copy()
    for species in pbesol_energy_eV:
        pbesol_energy_J[species] = pbesol_energy_eV[species]*eV2J*constants.N_A
    ### thermodynamic data for periodic solids ###
    thermo_file=dict()
    U_func = dict()

    phonon_species = ('czts','Cu','Zn','Sn','alpha_S')
    for species in phonon_species:
        thermo_file[species] = 'phonopy_output/' + species + '.dat'
        U_func[species] = get_potential_aims(thermo_file[species],'U')

    ### thermodynamic data for non-periodic components ###
    gas_species = ('S8')
    for species in gas_species:

        print 0

############## Enthalpy for each species #####################
    V = 1E-30 # Need to import correct cell volume from geometry file
    H = dict()
    for species in phonon_species:
        H[species] = pbesol_energy_eV[species] + U_func[species](T)*eV2J*constants.N_A + P*V

######### Calculate some formation energies! #################
    formation_alpha_eV = (pbesol_energy_eV['czts']/fu_cell['czts'] -
                          (2*pbesol_energy_eV['Cu']/fu_cell['Cu'] +
                           pbesol_energy_eV['Zn']/fu_cell['Zn'] +
                           pbesol_energy_eV['Sn']/fu_cell['Sn'] +
                           4*pbesol_energy_eV['alpha_S']/fu_cell['alpha_S']
                           ))
    print "Formation energy 2 Cu + Zn + Sn + 4 S(alpha) -> CZTS: {0:4.2f} eV".format(formation_alpha_eV)
                    
    return 0

if __name__=="__main__":
    main()
