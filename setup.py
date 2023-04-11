# -*- coding: utf-8 -*-
from setuptools import setup

packages = ['scripts']

package_data = {'': ['*']}

modules = ['detettore']

install_requires = \
['biopython==1.78',
 'numba==0.56.4',
 'numpy==1.21.5',
 'pysam==0.15.3',
 'scipy==1.7.3',
 'minimap2==2.24']

entry_points = \
{'console_scripts': ['bamstats = scripts.bamstats:main',
                     'combinevcf = scripts.combineVCFs:main',
                     'detettore = detettore:main']}

setup_kwargs = {
    'name': 'detettore',
    'version': '2.0.4',
    'description': 'A tool to detect transposable element polymorphisms',
    'long_description': None,
    'author': 'cstritt',
    'author_email': 'crstp.strt@gmail.com',
    'maintainer': 'cstritt',
    'maintainer_email': 'crstp.strt@gmail.com',
    'url': 'https://github.com/cstritt/detettore',
    'packages': packages,
    'package_data': package_data,
    'py_modules': modules,
    'install_requires': install_requires,
    'entry_points': entry_points,
    'python_requires': '==3.7',
}


setup(**setup_kwargs)
