import numpy as np
import matplotlib.pyplot as plt
from materials import Cu, Zn, Sn, alpha_S, CZTS, S8

T = np.linspace(0,1000,40)
P = 1E5

subplot=321

plt.figure(1)

for material in Cu, Zn, Sn, alpha_S, CZTS, S8:
    plt.subplot(subplot)
    plt.plot(T,material.mu_kJ(T,P))
    plt.xlabel('Temperature / T')
    plt.ylabel('Chemical potential / kJ mol-1 ')
    plt.title(material.name)
    subplot+=1

plt.tight_layout() # Prevent overlapping plot labels
plt.show()    


