import numpy as np
from materials import Cu, Zn, Sn, alpha_S, CZTS
from matplotlib import pyplot as plt
from matplotlib2tikz import save as tikz_save

T = np.linspace(0,1200,100)

for material in Cu, Zn, Sn, alpha_S:
    plt.plot(T,material.Cv_kB(T), label=material.name)

plt.plot(T,CZTS.Cv_kB(T)/8, label=CZTS.name)

plt.axis([0,1200,0,3])
plt.xlabel('Temperature (K)')
plt.ylabel('$C_v$ ($k_B$ atom$^{-1}$)')
plt.legend()
plt.savefig('Cv.pdf')
tikz_save("Cv.tex",figureheight="8cm",figurewidth="10cm")
#plt.show()
