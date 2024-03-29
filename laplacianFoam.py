#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: Bruno Ramoa
@affiliation: Institute for Polymers and Composites, University of Minho, Portugal
"""

import numpy as np
import sympy as sym
from sympy import cos

import pyMMSFoam as mms
from pyMMSFoam import x,y,z,t


# LaplacianFoam
T = 150*(cos(x*x + y*y + 0.1*t) + 1.5) # Unsteady 2D
# T = 150*(cos(x*x + y*y ) + 1.5) # Steady 2D
# T = 150*(cos(x*x + 0.1*t) + 1.5) # Unsteady 1D

DT = 1e-3

S_T = mms.ddt(T) - mms.laplacian(T, DT)


# # Uncomment what is needed

# # Generate fvOptions
# mms.generateFvOptions(S_T, "sourceTerm", "T")

# # Generate boundary conditions
# mms.generateDirichletBoundaries(T, "T")
# mms.generateNeumannBoundaries(T, "T")

# # Generate functionObjects
# mms.generateFunctionObject(T, "T")
                       
