#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: pc
"""

import sympy as sym
from .symbols import x,y,z

def checkInput(scalarOrVector):
    
     # Check if input is a int, float or symbolic expression
     checkIfScalar = isinstance(scalarOrVector, (int,float,sym.Expr))
     
     # Check if input is a vector (In this script vectors are [3,1] matrices)
     checkIfVectorOrTensor = isinstance(scalarOrVector, sym.matrices.dense.MutableDenseMatrix)
     
     
     checkIfVector = False
     checkIfTensor = False
     
     # Checks if the input is correctly formated
     if(checkIfScalar):
         checkScalar(scalarOrVector)
     elif(checkIfVectorOrTensor):
         if (scalarOrVector.shape == (3, 1)):
             checkVector(scalarOrVector)
             checkIfVector = True
         elif(scalarOrVector.shape == (3, 3)):
             checkTensor(scalarOrVector)
             checkIfTensor = True
         else:
             raise ValueError('There is a problem with the vector/tensor input')
     else:
         raise ValueError('Data structure not known')
         
     return [checkIfScalar, checkIfVector, checkIfTensor]
     

def checkScalar(S):
    
    # If S is a symbolic expression, check if x, y or z are part of the expression
    if (isinstance (S, sym.Expr)):
        try:
            float(S)
        except:
            assert( 
                    any( symbolCheck in S.free_symbols for symbolCheck in [x,y,z]) ), \
                        "Symbolic expressions must be defined in cartesian space with x,y and z \n " 

                  
def checkVector(V):  
 
    for i in V:
        try:
            float(i)
        except:
            assert( any( symbolCheck in V.free_symbols for symbolCheck in [x,y,z]) ), \
                    "Symbolic expressions must be defined in cartesian space with x,y and z" 

def checkTensor(T):   
                
    # Check if the entries of the tensor are symbolic expressions or numbers
    # Tries to convert each entry to a flot.
    # If it fails it assumes it is a symbolic expression and checks if x, y or z is present

    for i in T:
        try:
            float(i)
        except:
            assert( any( symbolCheck in T.free_symbols for symbolCheck in [x,y,z]) ), \
                "Symbolic expressions must be defined in cartesian space with x,y and z" 
                   