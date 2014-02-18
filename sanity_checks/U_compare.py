import numpy as np
from materials import Cu, Zn, Sn, alpha_S, CZTS_kesterite, CZTS_stannite
from materials import Cu2S_low, ZnS_zincblende, ZnS_wurtzite, SnS2, SnS_pcma
from matplotlib import pyplot as plt

T = np.linspace(0,1200,100)

for material in (Cu, Zn, Sn, alpha_S, CZTS_kesterite, CZTS_stannite,
                 Cu2S_low, ZnS_zincblende, ZnS_wurtzite, SnS2, SnS_pcma):
    plt.plot(T,(material.U_kJ(T)-material.U_kJ(0))/material.N,
             label=material.name)

#plt.axis([0,1200,0,3])
plt.xlabel('Temperature (K)')
plt.ylabel('$\Delta U$ (kJ atom mol$^{-1}$)')
plt.legend()
plt.savefig('plots/DeltaU.png')
plt.show()
