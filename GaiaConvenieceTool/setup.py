
from setuptools import setup, find_packages

setup(
    name='gaia',
    version='0.0.1',
    #url='https://github.com/mypackage.git',
    author='Ekaterina Ilin',
    author_email='eilin@aip.de',
    description='Some mini tools for gaia that I use.',
    packages=find_packages(),
    install_requires=['numpy', 'pandas'],
)