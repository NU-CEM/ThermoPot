import numpy as np

def find_stable_materials(potentials):

    assert (
        len(set([array.shape for array in potentials])) == 1
    ), "potential arrays must have the same dimension"

    minimum_potential = potentials[0]
    for i,potential in enumerate(potentials):
        minimum_potential = np.minimum(minimum_potential,potentials[i+1])
        if i + 2 == len(potentials):
            break

    for i,potential in enumerate(potentials):
        minimum_potential[potential == minimum_potential] = i

    return minimum_potential

