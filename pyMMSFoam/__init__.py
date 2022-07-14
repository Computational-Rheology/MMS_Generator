#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: pc
"""
from .symbols import x,y,z,t
from .checks import checkInput
from .differentialOperators import ddt, grad, div, laplacian
from .generateCCode import generateCcode
from .printUtilities import generateDirichletBoundaries, generateNeumannBoundaries, generateFunctionObject, generateFvOptions
from .plotting import plot
