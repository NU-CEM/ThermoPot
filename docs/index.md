ThermoPot is used for ab-initio thermodynamic modelling of material formation and decomposition.

- 📚 The documentation is [here](https://NU-CEM.github.io/ThermoPot). 
- 🔄 This code is made available under the GNU General Public Licence (GPL) v3. See the LICENSE file for the full text.
- 🖊 If you use this package for your research please [cite accordingly](https://github.com/NU-CEM/ThermoPot/blob/main/citation.cff).

This work adapts and extends a previous repository developed by [Adam Jackson](https://orcid.org/0000-0001-5272-6530): [Thermodynamic model of CZTS](http://dx.doi.org/10.5281/zenodo.57130). 

## Features

ThermoPot calculates temperature and pressure dependent thermodynamic potentials using first-principles data. Thermopot can:

- calculate the internal energy, enthalpy and Gibbs free energy at a given temperature and pressure
- calculate the change in energy/enthalpy for a given reaction
- work with solid or gases
- plot potentials as a function of T and P
- predict thermodynamic stability as a function of T and P
- Parse DFT data from an FHI-aims output file and lattice dynamics data from Phonopy output
- Parse experimental data from the NIST-JANAF thermochemical tables.

## Supported software

ThermoPot is compatible with a range of materials modelling packages as DFT energies and thermal properties data can be provided by the user. ThermoPot also supports parsing of FHI-aims output files, phonopy output and NIST-JANAF thermochemical data tables.

## Related software

There are other codes that can calculate thermodynamic potentials and phase diagrams. The best code depends on your use case:
- [Surfinpy](https://github.com/symmy596/SurfinPy) generates phase diagrams for bulk and surface systems. 
- [pycalphad](https://pycalphad.org/docs/latest/) uses the CALPHAD method to generate phase diagram and predict thermodynamic properties.
- [thermo](https://thermo.readthedocs.io/cubic_equations_of_state.html) uses cubic equations of state to generate phase diagrams for liquids and gases.
