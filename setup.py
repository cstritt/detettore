#!/usr/bin/env python

try:
    from setuptools import setup, find_packages
except ImportError:
    print '\nInstall first the Python setuptools module:\n\
    $ pip install setuptools, or with Anaconda:\n\
    $ conda install setuptools'

setup(name='detettore',
      version='1.1',
      description='A program to detect and characterize TE polymorphisms',
      author='Christoph Stritt',
      author_email='crstp.strt@gmail.com',
      url='https://github.com/cstritt/detettore',
      license='GPL',

      packages=find_packages(),
      
      scripts=['detettore.py',
               'filter.py',
               'variantcaller.py'],
      
               
      python_requires='<3',
               
      install_requires=[
              'pysam',
              'statistics',
              'biopython',
              'joblib',
              'scipy < 1.3',
              'setuptools'])

