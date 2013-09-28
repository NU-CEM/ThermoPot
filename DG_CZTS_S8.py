from materials import CZTS, Cu, Zn, Sn, S8
import matplotlib.pyplot as plt
import numpy as np

T = np.linspace(100,1000,10)    # K
P = np.array( np.logspace(3,7,10),ndmin=2).transpose() # Pa

D_mu = CZTS.mu_kJ(T,P) - (2*Cu.mu_kJ(T,P) +
                          Zn.mu_kJ(T,P) +
                          Sn.mu_kJ(T,P) +
                          0.5*S8.mu_kJ(T,P)
                          )


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
a = plt.contour(T,P.flatten(),D_mu,10, linewidths = 0.5, colors = 'k')
plt.pcolormesh(T,P.flatten(),D_mu,cmap=plt.get_cmap('BuPu'),vmin=-400, vmax=-220)
colours = plt.colorbar()
plt.xlabel('Temperature / K')
plt.ylabel('Pressure / Pa')
colours.set_label('$\Delta\mu_f / $kJ mol$^{-1}$')
ax.set_yscale('log')
plt.clabel(a)
plt.show()


