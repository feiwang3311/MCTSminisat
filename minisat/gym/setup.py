#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


example_module = Extension('_GymSolver',
                           sources=['GymSolver_wrap.cxx', 'GymSolver.cc'],
                           )

setup (name = 'GymSolver',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [example_module],
       py_modules = ["GymSolver"],
       )
