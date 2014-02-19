CZTS thermodynamic modelling
============================

Research data and calculations for ab initio thermodynamic modelling of
the formation of CZTS

Provided as supplementary information for publication, this is also an active project currently hosted at
[http://github.com/WMD-Bath/CZTS-model](http://github.com/WMD-Bath/CZTS-model).

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