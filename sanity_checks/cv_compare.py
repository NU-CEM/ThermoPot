import numpy as np
from materials import Cu, Zn, Sn, alpha_S, CZTS_kesterite, CZTS_stannite
from materials import Cu2S_low, ZnS_zincblende, ZnS_wurtzite, SnS2, SnS_pcma
from matplotlib import pyplot as plt

T = np.linspace(0,1200,100)

for material in (Cu, Zn, Sn, alpha_S, CZTS_kesterite, CZTS_stannite,
                 Cu2S_low, ZnS_zincblende, ZnS_wurtzite, SnS2, SnS_pcma):
    plt.plot(T,material.Cv_kB(T)/material.N, label=material.name)

plt.axis([0,1200,0,3])
plt.xlabel('Temperature (K)')
plt.ylabel('$C_v$ ($k_B$ atom$^{-1}$)')
plt.legend()
plt.savefig('plots/Cv.png')
plt.show()
