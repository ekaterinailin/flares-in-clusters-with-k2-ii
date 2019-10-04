# Flares in Clusters with K2. II. 
_Ilin+2019 (in prep.)_

The project is structured as follows:

1. *Membership matching* from different sources
2. *Stellar parameters* and quiescent luminosities from color-temperature relations and model spectra.
3. *Flare finding* and injection-recovery characterization of flare candidates with  with [*AltaiPony*](https://github.com/ekaterinailin/AltaiPony).
4. *Analysis* of flaring activity with respect to stellar mass and age.

Each part is enclosed in one folder that contains

- a README with instructions for those who wish to reproduce the results
- Jupyter notebooks that are introduced in the README file
- modules and scripts (Python 3.5)
- the resulting tables and plots
- relevant ancillary and raw data
- plots that are included in Ilin+2019 (in prep.)

NOTE: Some folders do not yet exist but will be added soon.

## Installation instructions

Clone this repository, create a virtual environment, activate it, then run the installation file, like so:

```
git clone https://github.com/ekaterinailin/flares-in-clusters-with-k2-ii.git
cd flares-in-clusters-with-k2-ii

python3 -m pip install --user virtualenv
python3 -m venv flaresinclustersii
source flaresinclustersii/bin/activate

bash installation.sh
```

Feel free to roam through the project, try the code, and browse the tables.

Questions, remarks, ideas for improvement? Open an issue in this repository or send an email to [Ekaterina Ilin](eilin@aip.de).
