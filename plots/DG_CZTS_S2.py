################################################################################
#  Copyright Adam J. Jackson (2014)                                            #
#                                                                              #
#   This program is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   This program is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
################################################################################

def main():
    from materials import CZTS, Cu, Zn, Sn, S2
    import numpy as np

    from DG_CZTS_S8 import plot_potential

    T = np.linspace(100,1500,100)    # K
    P = np.array( np.logspace(1,7,100),ndmin=2).transpose() # Pa
    
    D_mu = CZTS.mu_kJ(T,P) - (2*Cu.mu_kJ(T,P) +
                                    Zn.mu_kJ(T,P) +
                                    Sn.mu_kJ(T,P) +
                                    2*S2.mu_kJ(T,P)
        )
    
    D_mu_label = '$\Delta G_f$ / kJ mol$^{-1}$'
    scale_range = [-650,50]
    
    plot_potential(T,P,D_mu,D_mu_label,scale_range, filename='plots/DG_CZTS_S2.png')

if __name__ == "__main__":
    main()
