import numpy as np
from materials import Cu, Zn, Sn, alpha_S, CZTS_kesterite, CZTS_stannite
from matplotlib import pyplot as plt

T = np.linspace(0,1200,100)

for material in Cu, Zn, Sn, alpha_S:
    plt.plot(T,material.Cv_kB(T), label=material.name)

for material in CZTS_kesterite, CZTS_stannite:
    plt.plot(T,material.Cv_kB(T)/8, label=material.name)

plt.axis([0,1200,0,3])
plt.xlabel('Temperature (K)')
plt.ylabel('$C_v$ ($k_B$ atom$^{-1}$)')
plt.legend()
plt.savefig('plots/Cv.png')
plt.show()
