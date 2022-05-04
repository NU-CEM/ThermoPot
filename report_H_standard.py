from materials import BaZrS3, Ba2ZrS4, Ba3Zr2S7_RP, Ba3Zr2S7_Needle 
from materials import Ba, Zr, S, S2, S8
from materials import BaS, BaS2, BaS3, ZrS, ZrS2, ZrS3

T = 298.15  # K
P = 1E5     # Pa

#### Summary of possible plots

# The process used to find structures (funnel)
# Crystal structures of ternary structures, with charge densities alongside (lone pair formation)
# Phonon dispersion for BaZrS3
# Heats of formation for perovskite from solid elemental forms, and from gas sulfur (S2,S8)
# T-P diagram for perovskite from the stable binaries (ball milling)
# T-P diagram for perovskite from the stable binaries in sulfur gas rich environment (thin film annealing, nanoparticle)
# T-P diagram for: binaries to sulfur excess binaries in a sulfur rich environment (link to Scragg work)
# Gibbs vs T for: perovskite from stable binaries, binaries from sulfur excess binaries, the other ternary phases (all in sulfur rich)

## Heats of formation

#### BaZrS3 <---> Ba + Zr + 3S

DH_f_BaZrS3_eV_pbesol = BaZrS3.H_eV(T,P) - (Ba.H_eV(T,P) + Zr.H_eV(T,P)
                            + 3.* S.H_eV(T,P))
DH_f_BaZrS3_kJ_pbesol = BaZrS3.H_kJ(T,P) - (Ba.H_kJ(T,P) + Zr.H_kJ(T,P)
                            + 3.* S.H_kJ(T,P))

DE_f_BaZrS3_eV_pbesol = BaZrS3.pbesol_energy_eV/BaZrS3.fu_cell - (
          Ba.pbesol_energy_eV/Ba.fu_cell +
          Zr.pbesol_energy_eV/Zr.fu_cell +
          3.*S.pbesol_energy_eV/S.fu_cell
)

DH_f_BaZrS3_eV_hse06 = BaZrS3.H_eV(T,P,xc='hse06') - (Ba.H_eV(T,P,xc='hse06') + Zr.H_eV(T,P,xc='hse06')
                            + 3.* S.H_eV(T,P,xc='hse06'))
DH_f_BaZrS3_kJ_hse06 = BaZrS3.H_kJ(T,P,xc='hse06') - (Ba.H_kJ(T,P,xc='hse06') + Zr.H_kJ(T,P,xc='hse06')
                            + 3.* S.H_kJ(T,P,xc='hse06'))

DE_f_BaZrS3_eV_hse06 = BaZrS3.hse06_energy_eV/BaZrS3.fu_cell - (
          Ba.hse06_energy_eV/Ba.fu_cell +
          Zr.hse06_energy_eV/Zr.fu_cell +
          3.*S.hse06_energy_eV/S.fu_cell
)

print('###### PBEsol values #######')

print ('Formation enthalpy of kesterite CZTS: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_BaZrS3_eV_pbesol)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_BaZrS3_kJ_pbesol)
       )

print ('Ground-state formation energy:        ' +
       '{0:3.2f} eV / formula unit'.format(DE_f_BaZrS3_eV_pbesol)
       )

print ('Standard enthalpy - ground-state energy: ' +
       '{0:f} eV / formula unit'.format(DH_f_BaZrS3_eV_pbesol-DE_f_BaZrS3_eV_pbesol)
       )

print('###### HSE06 values #######')

print ('Formation enthalpy of kesterite CZTS: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_CZTS_eV_hse06)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_CZTS_kJ_hse06)
       )

print ('Ground-state formation energy:        ' +
       '{0:3.2f} eV / formula unit'.format(DE_f_CZTS_eV_hse06)
       )

print ('Standard enthalpy - ground-state energy: ' +
       '{0:f} eV / formula unit'.format(DH_f_CZTS_eV-DE_f_CZTS_eV_hse06)
       )

#### BaZrS3 <---> Ba + Zr + (3/2)*S2

# TODO

#### BaZrS3 <---> Ba + Zr + (3/8)*S8

# TODO

#### Ba2ZrS4 <--> 2Ba + Zr + 4*S

# TODO

#### Ba2ZrS4 <--> 2Ba + Zr + 2*S2

# TODO

#### Ba2ZrS4 <--> 2Ba + Zr + (1/2)*S8

# TODO

#### Ba3Zr2S7 <--> 2Ba + Zr + 4*S

# TODO

#### Ba3Zr2S7 <--> 2Ba + Zr + 2*S2

# TODO

#### Ba3Zr2S7 <--> 2Ba + Zr + (1/2)*S8

# TODO

## Heats of formation

#### BaZrS3 <---> Ba + Zr + 3S_s

"""Binary phases"""

#DH_f_Cu2S_eV = Cu2S.H_eV(T,P) - (
               2.*Cu.H_eV(T,P) + alpha_S.H_eV(T,P)
               )           

#DH_f_Cu2S_kJ = Cu2S.H_kJ(T,P) - (
               2.*Cu.H_kJ(T,P) + alpha_S.H_kJ(T,P)
               )           

#DH_f_SnS_eV = SnS.H_eV(T,P) - (
    Sn.H_eV(T,P) + alpha_S.H_eV(T,P)
               )

#DH_f_SnS_kJ = SnS.H_kJ(T,P) - (
    Sn.H_kJ(T,P) + alpha_S.H_kJ(T,P)
               )

#DH_f_SnS2_eV = SnS2.H_eV(T,P) - (
    Sn.H_eV(T,P) + 2.*alpha_S.H_eV(T,P)
               )

#DH_f_SnS2_kJ = SnS2.H_kJ(T,P) - (
    Sn.H_kJ(T,P) + 2.*alpha_S.H_kJ(T,P)
               )

#DH_f_ZnS_eV = ZnS.H_eV(T,P) - (
    Zn.H_eV(T,P) + alpha_S.H_eV(T,P)
               )

#DH_f_ZnS_kJ = ZnS.H_kJ(T,P) - (
    Zn.H_kJ(T,P) + alpha_S.H_kJ(T,P)
               )


print ('Formation enthalpy of low Cu2S: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_Cu2S_eV)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_Cu2S_kJ)
       )

print ('Formation enthalpy of SnS: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_SnS_eV)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_SnS_kJ)
       )

print ('Formation enthalpy of SnS2: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_SnS2_eV)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_SnS2_kJ)
       )


print ('Formation enthalpy of zinc blende: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_ZnS_eV)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_ZnS_kJ)
       )
