CZTS thermodynamic modelling
============================

Research data and calculations for ab initio thermodynamic modelling of
the formation and decomposition of Cu<sub>2</sub>ZnSnS<sub>4</sub> (CZTS).

Provided as supplementary information for publication, this is also an active project currently hosted at
[http://github.com/WMD-Bath/CZTS-model](http://github.com/WMD-Bath/CZTS-model).
The v1.0 tag indicates the version supplied with the project's first publication.

(c) Adam Jackson 2014
This code is made available under the GNU General Public Licence (GPL) v3.
See the LICENSE file for the full text.

Contents
--------

* **materials.py** Core python library containing:
  * Classes for material results, implementing key thermodynamic functions
  * Objects containing results data for materials by name:
    * DFT total energies
    * Basic structure parameters
    * Filenames of supporting data
    * Methods for calculation of T- and P-dependent thermodynamic potentials

* **interpolate_thermal_property.py** Interpolation functions for tabulated data, using Scipy.

* **nist_janaf/** Folder containing thermodynamic data for S2 and S8 gases from the literature.

* **phonopy_output/** Numerically tabulated thermodynamic properties from [Phonopy](http://phonopy.sourceforge.net) runs.

* **plots/** Plotting programs, primarily for computing free energy surfaces.
  The main plotting routine is contained in **plots/DG_CZTS_S8.py** and this function is imported
  as needed by other free energy surface plotting programs.
  Due to the structure of Python libraries, these functions need to be called from the parent folder, e.g.
  `python plots/DG_CZTS_S2.py`.

* **report_H_standard.py** Calculate and print key formation energies.

* **sanity\_checks** contains several functions which compare thermodynamic properties over
  a temperature range, rescaled for direct comparison. This is intended to highlight basic
  issues in calculations. The files need to be run from the parent directory, i.e.:

      python sanity_checks/mu_compare.py

* **jscragg_2011.csv** Data file containing kinetic model stability boundary
  data from [Scragg et al. (2011)](http://dx.doi.org/10.1021/cm202379s). This
  data is used in **plots/DG_CZTS_SnS_Scragg.py** to reproduce Fig. 7 of [Jackson and Walsh (2014)](http://dx.doi.org/10.1039/c4ta00892h).
