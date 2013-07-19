"""Calculate enthalpies of formation for CZTS over a range of temperatures at standard pressure."""
from scipy import constants # Library provides 2010 CODATA values for physical constants
from interpolate_phonopy import get_potential_aims

def main():
################ Set conditions ###############
    T = 298.15 # K
    P = 1e5    # Pa
############# Data and parameters from DFT calcs ################

    pbesol_energy_eV=dict(czts=-0.706480597450521e06,
                          Cu=-0.180838109862865e06,
                          Zn=-0.981595848099161e05, 
                          Sn=-0.681162478362528e06, 
                          alpha_S=-0.347575504588933e06,
                          S8=-0.868936310037924e05)
    
    fu_cell=dict(czts=2,
                 Cu=4,
                 Zn=2,
                 Sn=4,
                 alpha_S=32,
                 S8=1)
    
    thermo_file=dict()
    U_func = dict()

    phonon_species = ('czts','Cu','Zn','Sn','alpha_S')
    for species in phonon_species:
        thermo_file[species] = 'phonopy_output/' + species + '.dat'
        U_func[species] = get_potential_aims(thermo_file[species],'U')

    ### Copy energies and convert to J: ###
    eV2J = constants.physical_constants['electron volt-joule relationship'][0]
    pbesol_energy_J = pbesol_energy_eV.copy()
    for species in pbesol_energy_eV:
        pbesol_energy_J[species] = pbesol_energy_eV[species]*eV2J*constants.N_A

############## Enthalpy for each species #####################
        V = 1E-30 # Need to import correct cell volume from geometry file
    H = dict()
    for species in phonon_species:
        H[species] = pbesol_energy_eV[species] + U_func[species](T)*eV2J*constants.N_A + P*V
        print H

######### Calculate some formation energies! #################
    formation_alpha_eV = (pbesol_energy_eV['czts']/fu_cell['czts'] -
                          (2*pbesol_energy_eV['Cu']/fu_cell['Cu'] +
                           pbesol_energy_eV['Zn']/fu_cell['Zn'] +
                           pbesol_energy_eV['Sn']/fu_cell['Sn'] +
                           4*pbesol_energy_eV['alpha_S']/fu_cell['alpha_S']
                           ))
    print "Formation energy 2 Cu + Zn + Sn + 4 S(alpha) -> CZTS: {0:4.2f} eV".format(formation_alpha_eV)

    formation_alpha_J = (pbesol_energy_J['czts']/fu_cell['czts'] -
                         (2*pbesol_energy_J['Cu']/fu_cell['Cu'] +
                          pbesol_energy_J['Zn']/fu_cell['Zn'] +
                          pbesol_energy_J['Sn']/fu_cell['Sn'] +
                          4*pbesol_energy_J['alpha_S']/fu_cell['alpha_S']
                          ))
    print "Formation energy 2 Cu + Zn + Sn + 4 S(alpha) -> CZTS: {0:4.2f} kJ mol-1".format(formation_alpha_J/1000)

    formation_S8_eV = (pbesol_energy_eV['czts']/fu_cell['czts'] -
                       (2*pbesol_energy_eV['Cu']/fu_cell['Cu'] +
                        pbesol_energy_eV['Zn']/fu_cell['Zn'] +
                        pbesol_energy_eV['Sn']/fu_cell['Sn'] +
                        0.5*pbesol_energy_eV['S8']/fu_cell['S8']
                        ))
    print "Formation energy 2 Cu + Zn + Sn + 0.5 S8 -> CZTS: {0:4.2f} eV".format(formation_S8_eV)

    formation_S8_J = (pbesol_energy_J['czts']/fu_cell['czts'] -
                         (2*pbesol_energy_J['Cu']/fu_cell['Cu'] +
                          pbesol_energy_J['Zn']/fu_cell['Zn'] +
                          pbesol_energy_J['Sn']/fu_cell['Sn'] +
                          0.5*pbesol_energy_J['S8']/fu_cell['S8']
                          ))
    print "Formation energy 2 Cu + Zn + Sn + 0.5 S8 -> CZTS: {0:4.2f} kJ mol-1".format(formation_S8_J/1000)
                    
    return 0

if __name__=="__main__":
    main()
