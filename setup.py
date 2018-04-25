#!/usr/bin/env python3

from setuptools import setup

setup(name='genedescriptions',
      version='0.1',
      description='gene descriptions package',
      url='',
      author='Valerio Arnaboldi',
      author_email='valearna@caltech.edu',
      packages=['genedescriptions'],
      install_requires=[
          'namedlist',
          'inflect',
          'PyYAML',
          'numpy'
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
