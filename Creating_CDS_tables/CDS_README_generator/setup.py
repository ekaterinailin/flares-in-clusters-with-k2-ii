from setuptools import setup, find_packages

setup(
    name='cds',
    version='0.0.1',
    #url='https://github.com/mypackage.git',
    author='not me',
    author_email='someone@else.de',
    description='CDS README tool.',
    packages=find_packages(),
    install_requires=['numpy','astropy'],
    package_data={'': ['bytebybyte.template', 'ReadMe.template']},
)
