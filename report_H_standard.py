#! /usr/bin/env python
from materials import CZTS, Cu, Zn, Sn, alpha_S, Cu2S_low as Cu2S, SnS2, ZnS_zincblende as ZnS

T = 298.15  # K
P = 1E5     # Pa

DH_f_CZTS_eV = CZTS.H_eV(T,P) - (2.* Cu.H_eV(T,P) + Zn.H_eV(T,P)
                            + Sn.H_eV(T,P) + 4.* alpha_S.H_eV(T,P))
DH_f_CZTS_kJ = CZTS.H_kJ(T,P) - (2.* Cu.H_kJ(T,P) + Zn.H_kJ(T,P)
                            + Sn.H_kJ(T,P) + 4.* alpha_S.H_kJ(T,P))

DE_f_CZTS_eV = CZTS.pbesol_energy_eV/CZTS.fu_cell - (
          2.*Cu.pbesol_energy_eV/Cu.fu_cell +
          Zn.pbesol_energy_eV/Zn.fu_cell +
          Sn.pbesol_energy_eV/Sn.fu_cell +
          4.*alpha_S.pbesol_energy_eV/alpha_S.fu_cell
)

print ('Formation enthalpy of kesterite CZTS: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_CZTS_eV)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_CZTS_kJ)
       )

print ('Ground-state formation energy:        ' +
       '{0:3.2f} eV / formula unit'.format(DE_f_CZTS_eV)
       )

print ('Standard enthalpy - ground-state energy: ' +
       '{0:f} eV / formula unit'.format(DH_f_CZTS_eV-DE_f_CZTS_eV)
       )


"""Binary phases"""

DH_f_Cu2S_eV = Cu2S.H_eV(T,P) - (
               2.*Cu.H_eV(T,P) + alpha_S.H_eV(T,P)
               )           

DH_f_Cu2S_kJ = Cu2S.H_kJ(T,P) - (
               2.*Cu.H_kJ(T,P) + alpha_S.H_kJ(T,P)
               )           

DH_f_SnS2_eV = SnS2.H_eV(T,P) - (
    Sn.H_eV(T,P) + 2.*alpha_S.H_eV(T,P)
               )

DH_f_SnS2_kJ = SnS2.H_kJ(T,P) - (
    Sn.H_kJ(T,P) + 2.*alpha_S.H_kJ(T,P)
               )

DH_f_ZnS_eV = ZnS.H_eV(T,P) - (
    Zn.H_eV(T,P) + alpha_S.H_eV(T,P)
               )

DH_f_ZnS_kJ = ZnS.H_kJ(T,P) - (
    Zn.H_kJ(T,P) + alpha_S.H_kJ(T,P)
               )


print ('Formation enthalpy of low Cu2S: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_Cu2S_eV)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_Cu2S_kJ)
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
