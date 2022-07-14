#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: pc
"""

import sympy as sym
from .checks import checkInput
from .symbols import x,y,z,t

# Differential operators
def ddx(scalar):
    return sym.diff(scalar, x)

def ddy(scalar):     
    return sym.diff(scalar, y)

def ddz(scalar):
    return sym.diff(scalar, z)


def ddt(scalarOrVector):
    
    check = checkInput(scalarOrVector) 
        
    if (check[0]):
        # Scalar
        timeRateOfChange=sym.diff(scalarOrVector, t)
        
    elif(check[1]):
        # Vector
        timeRateOfChange=sym.Matrix([ 
                                      sym.diff(scalarOrVector[0], t),
                                    + sym.diff(scalarOrVector[1], t),
                                    + sym.diff(scalarOrVector[2], t), 
                                    ])
    return timeRateOfChange  


def gradScalar(scalar):
    checkInput(scalar)

    scalarGradient = sym.Matrix([
                                 ddx(scalar),
                                 ddy(scalar),
                                 ddz(scalar)
                             ])
    return scalarGradient

def divVector(vector):
    checkInput(vector)

    divergence = ( 
                          ddx(vector[0])
                        + ddy(vector[1])
                        + ddz(vector[2])
                     )
    return divergence


def grad(scalarOrVector):
    
    check = checkInput(scalarOrVector)
        
    if (check[0]):
        # Computes the gradient of a scalar
        gradient = gradScalar(scalarOrVector)
        
    elif(check[1]):
        # Computes the gradient of a vector according to OF syntax
        gradient = sym.Matrix([])
        for i in scalarOrVector:
            gradient = gradient.row_join(grad(i))
        
    else:
        raise ValueError('Some problem in the gradient operator')
        
    return gradient
    


def div(vectorOrTensor):
    
     # Check if the vector is in the correct format
    check = checkInput(vectorOrTensor)    
        
    if(check[1]):
        # Computes the divergence of a vector
        divergenceOfVectorOrTensor = divVector(vectorOrTensor)
    elif(check[2]):
        # Computes the divergence of a tensor to OF syntax
        divergenceOfVectorOrTensor = sym.Matrix([])
        for i in range(3):
            vector = vectorOrTensor[:,i]
            divVec = divVector(vector) 
            divergenceOfVectorOrTensor = divergenceOfVectorOrTensor.col_join(sym.Matrix([divVec]))
    else:
        raise ValueError('Some problem in the divergence operator')
        
    return divergenceOfVectorOrTensor

def laplacian(S ,coef=1):    
    
    # Check if the scalar is in the correct format 
    check = checkInput(S)
        
    if (check[0]):
        laplcian = div(coef*grad(S))
    else:
        raise ValueError('Some problem in the laplacian operator')
    
    return laplcian

# def curl(v):    
    
#     # Check if the vector is in the correct format 
#     check = checkInput(v)
        
#     if (check[1]):
#         curl_ = sym.Matrix([
#                                 ddy(v[2])-ddz(v[1]),
#                                 ddz(v[0])-ddx(v[2]),
#                                 ddx(v[1])-ddy(v[0])
#                           ])
#     else:
#         raise ValueError('Some problem in the curl operator. Check If input is a vector')
    
#     return curl_


