#bash

# You should run this within a virtualenv, like so:
# python3.5
#python3 -m pip install --user virtualenv
#python3 -m venv flaresinclustersii
#source flaresinclustersii/bin/activate

#this now should run within the virtualenv

# get some pip installable packages
pip install dustmaps 
pip install bokeh

# get others that are not
mkdir SpecMatchEmp
cd SpecMatchEmp
git clone https://github.com/samuelyeewl/specmatch-emp
cd specmatch-emp
python setup.py install

cd ../..

#get AltaiPony 
mkdir AltaiPony
cd AltaiPony
git clone https://github.com/ekaterinailin/altaipony
cd altaipony
python setup.py install

cd ../..

# and k2sc
mkdir K2SC
cd K2SC
git clone https://github.com/ekaterinailin/k2sc.git
cd k2sc
python3 setup.py install

# install your own
cd ../../FlareAnalysisPipeline
python setup.py install
cd ../GaiaConvenienceTool
python setup.py install
cd ../Creating_CDS_tables/CDS_README_generator
python setup.py install
cd ../..
## to close the virtualenv again:
#deactive
