#bash

# python3.5
python3 -m pip install --user virtualenv
python3 -m venv flaresinclustersii
source flaresinclustersii/bin/activate

#this now should run within the virtualenv

# get some pip installable packages
pip install dustmaps 
pip install bokeh

# get others that are not
git clone https://github.com/samuelyeewl/specmatch-emp
python specmatch-emp/setup.py install

# install your own
python FlareAnalysisPipeline/setup.py install
python GaiaConvenienceTool/setup.py install
python Creating_CDS_tables/CDS_README_generator/setup.py install

## to close the virtualenv again:
#deactive
