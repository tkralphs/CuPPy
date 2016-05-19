#!/usr/bin/env python

from setuptools import setup
import setuptools

setup(name='coinor.cuppy',
      version='0.7.1',
      description='Cutting Plane Methods in Python (CuPPy)',
      long_description='''A collection of 'naive' implementations of basic
      cutting plane algorithms in Python. The collection contains a generator
      for Gomory Mixed Integer cuts and one for generating the most violated
      split cut using the method of Saxena and Balas.

      The idea is for these implementations to be as transparent as possible.
      I intend them mainly for educational use. They will most certainly not
      be effective in a real-world environment. Even on small examples, it is
      easy to run into numerical difficulties.

      The underlying solvers are Clp (https://projects.coin-or.org/Clp) and
      Cbc (https://projects.coin-or.org/Cbc), which is called via the Python
      bindings of CyLP (https://github.com/coin-or/CyLP). Through CyLP, one
      can easily add these as cut generators within Cbc.

      The cutting plane procedure can be visualized for 2D examples using the
      polyhedron2D class of GrUMPy (https://github.com/coin-or/GrUMPy).''',
      author='Ted Ralphs',
      author_email='ted@lehigh.edu',
      license='Eclipse Public License',
      url='https://github.com/tkralphs/CuPPy/',
      namespace_packages=['coinor'],
      packages=['coinor.cuppy','coinor.cuppy.examples','coinor'],
      package_dir = {'coinor': 'src'},
      install_requires=['coinor.grumpy>=0.7.0', 'pulp', 'cylp>=0.7.0']
     )
