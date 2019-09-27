from setuptools import setup, find_packages

setup(
    name='opencluster',
    version='0.0.1',
    #url='https://github.com/mypackage.git',
    author='Ekaterina Ilin',
    author_email='eilin@aip.de',
    description='Analysis for Flares in Clusters II study',
    packages=find_packages(),
    install_requires=['numpy', 'pandas'],
    package_data={'': ['static/*', 'tests/tesfiles/*', 'clusters/*']},
)
