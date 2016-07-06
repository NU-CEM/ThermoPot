CZTS thermodynamic modelling
============================

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.18732.svg)](http://dx.doi.org/10.5281/zenodo.18732)

Research data and calculations for ab initio thermodynamic modelling of
the formation and decomposition of Cu<sub>2</sub>ZnSnS<sub>4</sub> (CZTS).

This repository acts as supplementary information for a [2014 publication in *J. Mater. Chem. A*](http://dx.doi.org/10.1039/C4TA00892H), and is also an active project currently hosted at
[http://github.com/WMD-Group/CZTS-model](http://github.com/WMD-Group/CZTS-model).
The [releases](https://github.com/WMD-group/CZTS-model/releases) correspond to key publication points in the project:

* [report-confirmation](https://github.com/WMD-group/CZTS-model/releases/tag/report-confirmation) [ajjackson](https://github.com/ajjackson)'s 1st year PhD confirmation report. [Not public]
* [v1.0](https://github.com/WMD-group/CZTS-model/releases/tag/v1.0) Initial submission to *J. Mater. Chem. A*.
* [v1.2a](https://github.com/WMD-group/CZTS-model/releases/tag/v1.2a) Supporting data for publication in *J. Mater. Chem. A.*. Includes minor bug fixes and data for comparison with [another study](https://dx.doi.org/10.1021/cm202379s), with permission from Jonathan Scragg.
* [thesis-submission](https://github.com/WMD-group/CZTS-model/releases/tag/thesis-submission) Supporting data for ajjackson PhD thesis submission. [Initial submission is not public]




(c) Adam Jackson 2016
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

Example
-------

Each material is made available as an object in the module namespace,
and has methods which retrieve various properties. For example, to
generate a CSV file containing the chemical potentials of a selection of phases:

``` python
import csv
import numpy as np
from materials import CZTS_kesterite, Cu, beta_Sn, alpha_Sn
from materials import Cu2S, SnS2, SnS, Sn2S3, ZnS, S2, S8

T = np.arange(400, 1200, 50)

titles = []
data = []
for material in (CZTS_kesterite, Cu, beta_Sn,
                 alpha_Sn, Cu2S, SnS2,
                 SnS, Sn2S3, ZnS, S2, S8):

    titles.append(material.name)
    data.append(list(material.mu_kJ(T, 1E5)))

# Zip and list unpacking can be used together to transpose
# a matrix expressed as a list of lists
data = zip(*data)

with open('mu_data_Jmol.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerow(titles)
    writer.writerows(data)
```
