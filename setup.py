#!/usr/bin/env python
# encoding: utf-8
import sys
from distutils.command.build_py import build_py
from radexgrid import __version__ as version

if 'develop' in sys.argv:
    from setuptools import setup
else:
    from distutils.core import setup


with open('README.rst', 'r') as readme:
    long_description = readme.read()

setup(name='radexgrid',
      version=version,
      description='radexgrid - Python wrapper for RADEX',
      long_description=long_description,
      author=['Brian Svoboda'],
      author_email=['svobodb@email.arizona.edu'],
      url='https://github.com/autocorr/besl',
      packages=['radexgrid'],
      provides=['radexgrid'],
      requires=['numpy', 'pandas', 'astropy'])


