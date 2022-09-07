#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: pc
"""

# Script for MMS in OpenFOAM
import sympy as sym
from sympy import sin, cos, exp, pi, sqrt
import pyMMSFoam as mms
from pyMMSFoam import x,y,z,t

Re = 5
Lambda = (Re/2) - sqrt( (Re**2/4) + 4*pi**2 )
u = 1 - exp(Lambda*x)*cos(2*pi*y)
v = (Lambda/(2*pi))*exp(Lambda*x)*sin(2*pi*y)
w = 0
p = 0.5*(1-exp(2*Lambda*x))

U = sym.Matrix([u,v,w])

# Momentum balance equation
nu = 0.01

R = nu*(mms.grad(U) + mms.grad(U).T)

S = mms.div(U*U.T) - mms.div(R) + mms.grad(p)

# # Uncomment what is needed

# # Generate fvOptions
# mms.generateFvOptions(S, "momentumSource", "U")

# # Generate boundary conditions
# # Velocity
# mms.generateDirichletBoundaries(U, "U")
# mms.generateNeumannBoundaries(U, "U")

# # Pressure
# mms.generateDirichletBoundaries(p, "p")
# mms.generateNeumannBoundaries(p, "p")

# # Generate functionObjects
# mms.generateFunctionObject(U, "U")
# mms.generateFunctionObject(p, "p")