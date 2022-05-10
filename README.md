Ab-initio thermodynamic modelling package
============================

Software for ab-initio thermodynamic modelling of material formation and decomposition.

Code examples given for the BaZrS<sub>3</sub> system.

*Note: This repository has been built on one previously developed for a [thermodynamic model of CZTS](http://dx.doi.org/10.5281/zenodo.57130). 
It adapts and extends (with copyright permissions) code originally designed and written by Adam Jackson.*

This code is made available under the GNU General Public Licence (GPL) v3.
See the LICENSE file for the full text.

Contents
--------

* **materials.py** Core python library containing:
  * Classes for material results, implementing key thermodynamic functions
  * Objects containing results data for materials by name:
    * DFT total energies (PBEsol)
    * DFT total energies (HSE06)
    * Basic structure parameters
    * Filenames of supporting data
    * Methods for calculation of T- and P-dependent thermodynamic potentials

* **interpolate_thermal_property.py** Interpolation functions for tabulated data, using Scipy.

* **structures/** Folder containing structural models for materials.

* **nist_janaf/** Folder containing thermodynamic data for S2 and S8 gases from the literature.

* **phonopy_output/** Numerically tabulated thermodynamic properties from [Phonopy](http://phonopy.github.io/phonopy/) runs.

* **plots/** Plotting programs, primarily for computing free energy surfaces and Ellingham diagrames.
  The main plotting routine is contained in **(TODO)** and this function is imported
  as needed by other free energy surface plotting programs.
  Due to the structure of Python libraries, these functions need to be called from the parent folder, e.g.
  `python plots/DG_BaZrS3_S2.py`.

* **report_H_standard.py** Calculate and print key formation energies.


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
