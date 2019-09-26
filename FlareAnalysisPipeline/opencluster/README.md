# OpenCluster

This is the the pre- and post- analysis for the _Flaring activity in open clusters with K2_ project.

The modules are grouped by _Core_,  _Post-analysis_, _Ancillary_, _Helper_, and _Paper_.


## Core modules

### ffd.py

FFD class

### opencluster.py

Main OpenCluster class



## Post-analysis modules

### analysis.py

### results.py

Results class

### statistics.py


## Ancillary modules

### cmd.py

Interactive CMDs that allows us to click-and-go flag outliers in different color-magnitude spaces.

### linfit.py

### lum.py

### powerlawfit.py

## Helper packages and modules

### get_pip.py

This is just here for reference. Colab does not have pip3 installed. running get_pip.py fixes that.

### layout.py

Contains colors, markers, string representations useful for plotting.

### specmatchemp/

Has its own README with details. Contains the Emprical Spectral Library with model-independent stellar parameters covering $T_\mathrm{eff}=3000-7000$ K. We use it to fit spectra to our target to get a better luminosity.

### utils.py

## Paper modules

### paperplots.py

### papertables.py


## Tests

