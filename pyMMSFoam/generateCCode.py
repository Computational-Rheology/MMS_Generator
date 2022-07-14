#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: pc
"""


  
# Generate C code from MMS
import sympy as sym 
from sympy.codegen.rewriting import create_expand_pow_optimization
import re


def vectorSyntax (solutionName, isTensor=False):
    if(not(isTensor)):
        return  ["d"+solutionName+"_dx" ,  "d"+solutionName+"_dy" , "d"+solutionName+"_dz" ]
    else:
        return [solutionName + "_xx", solutionName + "_xy", solutionName + "_xz",
                solutionName + "_yx", solutionName + "_yy", solutionName + "_yz",
                solutionName + "_zx", solutionName + "_zy", solutionName + "_zz"]
    
def standardSyntax (solutionName, isTensor=False):
    if(not(isTensor)):
        return  [solutionName + "_1", solutionName + "_2" , solutionName + "_3"]
    else:
        return [solutionName + "_11", solutionName + "_12", solutionName + "_13",
                solutionName + "_21", solutionName + "_22", solutionName + "_23",
                solutionName + "_31", solutionName + "_32", solutionName + "_33"]
        
def getValue(isScalar, isVector):
    if (isScalar):
        return "0"
    elif(isVector):
        return "(0 0 0)"
    else:
        raise ValueError('variable in not defined correctly. \n Should be either scalar or vector')


def getTypeOfField(isScalar, isVector):
    if (isScalar):
        return "scalarField"
    elif(isVector):
        return "vectorField"
    else:
        raise ValueError('variable in not defined correctly.')


# Function to generate C code using common subexpression elimination for shortning the amount of code
def generateCcode(expression, solutionName, tempVarName='tmp', vectorNotation=False):
    
    # isScalar = isinstance(expression, sym.Expr)
    isVector = ( isinstance(expression, sym.matrices.dense.MutableDenseMatrix) and expression.shape == (3, 1) )
    isTensor = ( isinstance(expression, sym.matrices.dense.MutableDenseMatrix) and expression.shape == (3, 3) )
    

    if(vectorNotation):
        if(isVector):
            solutionName = vectorSyntax(solutionName)
        elif(isTensor):
            solutionName = vectorSyntax(solutionName, isTensor = True)
    else:
        if(isVector):
            solutionName = standardSyntax(solutionName)
        elif(isTensor):
            solutionName = standardSyntax(solutionName, isTensor = True)

    functions=['sin','cos','pow','exp','tan', 'sqrt', 'log']    

    code = sym.cse(expression, sym.numbered_symbols(tempVarName), optimizations='basic')
    expand_pow_if_needed = create_expand_pow_optimization(6)
    
    # Creates a list to store the results
    lines = []
    for tmp in code[0]:
                lines.append(sym.ccode(expand_pow_if_needed(tmp[1]), tmp[0]))

    for i, result in enumerate(code[1]):
        lines.append(sym.ccode(expand_pow_if_needed(result),  solutionName))
    
    if(isVector or isTensor):
        # gets the solution string 
        solution = lines[-1]
        # stips the solution from the list
        lines = lines[:-1]
        
        # strips the solution by new line caracter
        tmp_1 = solution.split("\n")
        
        # Adds solution back in the list for printing
        for i in tmp_1:
            lines.append(i)
    
    # Adds Foam:: to each function in the variable 'function' if in the string
    for i in range(len(lines)):   
        for j in range(len(functions)):
            lines[i] = re.sub(r'\b{}\b'.format(functions[j]), "Foam::"+functions[j], lines[i])
        lines[i] = "{0}const scalar " + lines[i] + " \n"

    return "".join(lines)    
