import numpy as np
from materials import Cu, Zn, Sn, alpha_S, CZTS_kesterite, CZTS_stannite
from materials import Cu2S_low, ZnS_zincblende, ZnS_wurtzite, SnS2, SnS_pcma
from materials import Cu2SnS3_mo1, Cu2SnS3_mo2
from matplotlib import pyplot as plt

from itertools import cycle
linestyles = ["-","--","-."]
linecycler = cycle(linestyles) # Set linecycler

T = np.linspace(0,1200,100)
P = 1E5                      # Pa

material_names=('Cu', 'Zn', 'Sn', r'$\alpha-S$', 'CZTS (kesterite)',
                'CZTS (stannite)', r'Low Cu$_2$S', 'ZnS (zincblende)',
                'ZnS (wurtzite)', r'SnS$_2$', 'SnS (pcma)',
                r'Cu$_2$SnS$_3$ (Mo-1)', r'Cu$_2$SnS$_3$ (Mo-2)')

fig = plt.figure(figsize=(7,5), dpi=200)

lines=[]
axis_list=[]

axis_list.append(plt.subplot(1,2,1))

for material in (Cu, Zn, Sn, alpha_S, CZTS_kesterite, CZTS_stannite,
                 Cu2S_low, ZnS_zincblende, ZnS_wurtzite, SnS2, SnS_pcma,
                 Cu2SnS3_mo1, Cu2SnS3_mo2):
    line, = plt.plot(T,material.Cv_kB(T)/material.N, next(linecycler), label=material.name)
    lines.append(line)

plt.axis([0,1200,0,3])
plt.xlabel('Temperature (K)')
plt.ylabel('$C_v$ ($k_B$ atom$^{-1}$)')

axis_list.append(plt.subplot(1,2,2))
#linecycler = cycle(lines) # Reset linecycler
for material in (Cu, Zn, Sn, alpha_S, CZTS_kesterite, CZTS_stannite,
                 Cu2S_low, ZnS_zincblende, ZnS_wurtzite, SnS2, SnS_pcma,
                 Cu2SnS3_mo1, Cu2SnS3_mo2):
    plt.plot(T,(material.mu_kJ(T,P)-material.mu_kJ(0,P))/material.N,
             next(linecycler), label=material.name)
plt.xlabel('Temperature (K)')
plt.ylabel('$\Delta \mu$ (kJ atom mol$^{-1}$)')

# Fix horizontal spacing
fig.tight_layout()
# Shrink axis heights by 30% on the bottom
for ax in axis_list:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.30,
                        box.width, box.height * 0.70])

plt.figlegend(lines, material_names,
              loc='upper center', fontsize='small',
              bbox_to_anchor=(0.5,0.25), ncol=3, fancybox=True,
              shadow=True
              )

plt.savefig('plots/Cv_mu_check.png')
plt.show()
