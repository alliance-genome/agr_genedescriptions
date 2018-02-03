#!/usr/bin/env python3

from setuptools import setup

setup(name='wb_genedescriptions',
      version='0.1',
      description='Generate gene descriptions',
      url='',
      author='Valerio Arnaboldi',
      author_email='valearna@caltech.edu',
      install_requires=[
          'namedlist'
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
