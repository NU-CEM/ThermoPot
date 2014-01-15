#! /usr/bin/env python
from materials import CZTS, Cu, Zn, Sn, alpha_S

T = 298.15  # K
P = 1E5     # Pa

DH_f_eV = CZTS.H_eV(T,P) - (2.* Cu.H_eV(T,P) + Zn.H_eV(T,P)
                            + Sn.H_eV(T,P) + 4.* alpha_S.H_eV(T,P))
DH_f_kJ = CZTS.H_kJ(T,P) - (2.* Cu.H_kJ(T,P) + Zn.H_kJ(T,P)
                            + Sn.H_kJ(T,P) + 4.* alpha_S.H_kJ(T,P))

DE_f_eV = CZTS.pbesol_energy_eV/CZTS.fu_cell - (
          2.*Cu.pbesol_energy_eV/Cu.fu_cell +
          Zn.pbesol_energy_eV/Zn.fu_cell +
          Sn.pbesol_energy_eV/Sn.fu_cell +
          4.*alpha_S.pbesol_energy_eV/alpha_S.fu_cell
)

print ('Formation enthalpy of kesterite CZTS: ' +
       '{0:3.2f} eV / formula unit'.format(DH_f_eV)
       )

print ('                                      ' +
       '{0:3.2f} kJ / mol'.format(DH_f_kJ)
       )

print ('Ground-state formation energy:        ' +
       '{0:3.2f} eV / formula unit'.format(DE_f_eV)
       )

print ('Standard enthalpy - ground-state energy: ' +
       '{0:f} eV / formula unit'.format(DH_f_eV-DE_f_eV)
       )
