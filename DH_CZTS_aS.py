"""Calculate enthalpies of formation for CZTS over a range of temperatures at standard pressure."""

def main():

    pbesol_energy_eV=dict(czts=-0.706480597450521e06,
                          Cu=-0.180838109862865e06,
                          Zn=-0.981595848099161e05, 
                          Sn=-0.681162478362528e06, 
                          alpha_S=-0.347575504588933e06)
    
    fu_cell=dict(czts=2,
                 Cu=4,
                 Zn=2,
                 Sn=4,
                 alpha_S=32)

    formation_eV = (pbesol_energy_eV['czts']/fu_cell['czts'] -
                    (2*pbesol_energy_eV['Cu']/fu_cell['Cu'] +
                     pbesol_energy_eV['Zn']/fu_cell['Zn'] +
                     pbesol_energy_eV['Sn']/fu_cell['Sn'] +
                     4*pbesol_energy_eV['alpha_S']/fu_cell['alpha_S']
                     ))
    DHf = formation_eV;

    return DHf

if __name__=="__main__":
    main()
