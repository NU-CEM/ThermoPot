""" Computed values for CZTS and related materials. """


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
                 alpha_S=32)

#### MADE-UP VOLUMES TO BE REPLACED!!! ###                 
volume=dict(czts=310.86645888987351,
                 Cu=45.38855878494433,
                 Zn=27.956078493840639,
                 Sn=108.8943337793184,
                 alpha_S=832.91786077871541
                 )


def volume_calc(filename):
    """Calculate unit cell volume in cubic angstroms from geometry.in file"""
    import numpy as np
    lattice_vectors = []
    with open(filename, 'r') as f:
        for line in f:
            if line.split()[0] == 'lattice_vector':
                lattice_vectors.append(line.split()[1:4])

    lattice_vectors = np.array(lattice_vectors).astype(float)
    volume = np.dot(lattice_vectors[0],np.cross(lattice_vectors[1],lattice_vectors[2]))

    return abs(volume)

    
