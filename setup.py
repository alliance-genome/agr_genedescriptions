#!/usr/bin/env python3

from setuptools import setup

setup(name='genedescriptions',
      version='2.0.7',
      description='gene descriptions package',
      url='',
      author='Valerio Arnaboldi',
      author_email='valearna@caltech.edu',
      packages=['genedescriptions'],
      install_requires=[
          'namedlist',
          'inflect',
          'PyYAML',
          'numpy',
          'urllib3',
          'ontobio'
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
