#!/usr/bin/env python

# this setup file does not work. See the ../../Makefile target python-wrap:

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


example_module = Extension('_GymSolver',
                           sources=['GymSolver_wrap.cc', 'GymSolver.cc'], include_dirs = ['/home/fei/progsyn/Reinforcement Learning/PolicyGradient/baselines/baselines/MCTS/minisat/'],
                           )

setup (name = 'GymSolver',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [example_module],
       py_modules = ["GymSolver"],
       )
