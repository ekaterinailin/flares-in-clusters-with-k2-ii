
from setuptools import setup, find_packages

setup(
    name='gaia',
    version='0.0.1',
    #url='https://github.com/mypackage.git',
    author='Ekaterina Ilin',
    author_email='eilin@aip.de',
    description='Some mini tools for gaia that I use.',
    package_data={'': ['readme.txt','table_u0_2D.txt','table_u0_g_col.txt','table_u0_g.txt']},
    packages=find_packages(),
    install_requires=['numpy', 'pandas'],
)
