#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='detettore',
      version='1.3',
      description='A program to detect and characterize TE polymorphisms',
      author='Christoph Stritt',
      author_email='crstp.strt@gmail.com',
      url='https://github.com/cstritt/detettore',
      license='GPL',

      packages=find_packages(),

      scripts=['detettore.py',
               'filter.py',
               'variantcaller.py'],

      python_requires='>=3.7',

      install_requires=[
              'pysam',
              'statistics',
              'biopython',
              'joblib',
              'scipy',
              'setuptools'])
