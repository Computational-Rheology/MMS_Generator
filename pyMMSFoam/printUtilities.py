#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: Bruno Ramoa
@affiliation: Institute for Polymers and Composites, University of Minho, Portugal
"""


# Checks
import sympy as sym 
from .generateCCode import generateCcode, getValue, getTypeOfField, standardSyntax, vectorSyntax
from .differentialOperators import grad
from .symbols import x,y,z,t

def getDimensions(expression):
    
    symbolsInExpression = expression.free_symbols
    
    tmp_dimensions =""
    
    x_definition = (
                        "{0}// Gets the x component of the current cell \n"
                        "{0}const scalar x = {1}{2}.x(); \n"
                    )
  
    y_definition =  (
                        "\n{0}//Gets the y component of the current cell \n"
                        "{0}const scalar y = {1}{2}.y(); \n "
                     )
 
    z_definition =(
                    "\n{0}// Gets the z component of the current cell \n"
                    "{0}const scalar z = {1}{2}.z(); \n"
                    )

        
    
    if(x in symbolsInExpression):
        tmp_dimensions = tmp_dimensions + x_definition
    if(y in symbolsInExpression):
        tmp_dimensions = tmp_dimensions + y_definition
    if(z in symbolsInExpression):
        tmp_dimensions = tmp_dimensions + z_definition
        
    return tmp_dimensions


def getTransient(expression, fvOptions = False):
    
    symbolsInExpression = expression.free_symbols
    
    tmp_transient =""
    if(t in symbolsInExpression):
        if(not(fvOptions)):
            tmp_transient = "{0}const scalar t = this->db().time().value(); \n"
        else:
            tmp_transient =( "\n{0}const Time& time = mesh().time(); \n"
                             "{0}// Gets the current time value \n"
                             "{0}const scalar t = time.value(); \n"
                             )
    return tmp_transient


def generateDirichletBoundaries(MMS, solutionName, tempVarName ='tmp'):
    
    isScalar = isinstance(MMS, sym.Expr)
    isVector = ( isinstance(MMS, sym.matrices.dense.MutableDenseMatrix) and MMS.shape == (3, 1) )

    value = getValue(isScalar, isVector)

    name = solutionName + "_dirichlet"
   
    transient = getTransient(MMS)
    
    dimensions = getDimensions(MMS)
    
    typeOfField = getTypeOfField(isScalar, isVector)

    solution = ""
    
    if(isScalar):
        solution = "field[faceI] = " + solutionName + " ;"

    elif (isVector):
        solNames = solutionName+"_1" +", " + solutionName +"_2" + ", " + solutionName+"_3"
        solution = (
                        "\n\t\t\tconst vector " + solutionName +"( "+ solNames + " ); \n"
                        "\t\t\tfield[faceI] = " + solutionName + "; \n"
                    )
    else:
        raise ValueError('variable in generateDirichletBoundaries is not defined correctly. \n Should be either scalar or vector')


    mainBody = (
                    "\"(patch1|patch2|patch3)\"                     \n"     
                    "{                                              \n"     
                    "\t// Dirichlet boundary                        \n"     
                    "\ttype        codedFixedValue;                 \n"      
                    "\tvalue       uniform %s;                      \n" 
                    "                                               \n" 
                    "\tname        %s;                              \n"  
                    "                                               \n" 
                    "\tcode                                         \n" 
                    "\t#{                                           \n" 
                    "                                               \n"
                    "\t\t// Gets current patch                      \n"
                    "\t\tconst fvPatch& boundaryPatch = patch();    \n"
                    "                                               \n"
                    "\t\t// Gets the patch face centres values      \n"
                    "\t\tconst vectorField& Cf = boundaryPatch.Cf(); \n"
                    "                                               \n"
                    "\t\t// Gets the current field                  \n"
                    "\t\t%s& field = *this;                         \n"
                    "                                               \n"
                    "\t\t// MMS                                     \n"
                    "%s                                             \n"
                    "\t\t// Loops over the patch                    \n"
                    "\t\tforAll(Cf, faceI)                          \n"
                    "\t\t{                                          \n"
                    "%s                                             \n"
                    "%s                                             \n"
                    "\t\t\t%s                                       \n"
                    "\t\t}                                          \n\n"
                    "\t#};                                          \n"
                    "}"
               )
    
    
    print(mainBody
                      %(
                          value,
                          name,
                          typeOfField,
                          transient.format("\t\t"),
                          dimensions.format("\t\t\t", "Cf","[faceI]"),
                          generateCcode(MMS, solutionName, tempVarName).format("\t\t\t"),
                          solution
                        )
         )
                 

def generateNeumannBoundaries(MMS, solutionName, tempVarName = 'tmp'):
    
    isScalar = isinstance(MMS, sym.Expr)
    
    isVector = ( isinstance(MMS, sym.matrices.dense.MutableDenseMatrix) and MMS.shape == (3, 1) )
    
    value = getValue(isScalar, isVector)
    
    name = solutionName + "_Neumann"
    
    transient = getTransient(MMS)
    
    dimensions = getDimensions(MMS)

    solNames = ""
    solution = ""
    
    if(isScalar):
        solNames =(  
                    "d" + solutionName + "_dx" +
                    ", " +
                    "d" + solutionName + "_dy" +
                    ", " +
                    "d" + solutionName + "_dz"
                  )
        solution = (
                "\n\t\t\tconst vector grad" + solutionName +" (" + solNames + ");   \n"
                "\n\t\t\tconst scalar normalGradient = " + "grad" + solutionName + " & nf[faceI] ; \n"
                "\t\t\tthis->refGrad()[faceI] = normalGradient; \n"
                "\t\t\tthis->valueFraction()[faceI] = scalar(0); \n"
                 )
        
    elif (isVector):
        solNames = [
                    solutionName+"1", solutionName+"2", solutionName+"3"
                    ]
        
        
        solution = (
                "\n\t\t\tconst scalar normal_1 = vector(" + ", ".join(vectorSyntax(solNames[0])) + ") & nf[faceI] ; \n"
                "\n\t\t\tconst scalar normal_2 = vector(" + ", ".join(vectorSyntax(solNames[1])) + ") & nf[faceI] ; \n"
                "\n\t\t\tconst scalar normal_3 = vector(" + ", ".join(vectorSyntax(solNames[2])) + ") & nf[faceI] ; \n"               
                
                "\n\t\t\tconst vector normalGradient (normal_1, normal_2, normal_3); \n"
                
                "\t\t\tthis->refGrad()[faceI] = normalGradient; \n"
                
                "\t\t\tthis->valueFraction()[faceI] = scalar(0); \n"
                 )
    else:
        raise ValueError('variable in generateNeumannBoundaries is not defined correctly. \n Should be either vector or tensor')
    
    code = []
    if (isScalar):              
        code.append(generateCcode(grad(MMS), solutionName, tempVarName, vectorNotation=True))
    elif(isVector):
        for i in range(len(solNames)):
            code.append(generateCcode(grad(MMS[i]), solNames[i], str("tmp" + str(i) + "_") , vectorNotation=True))
    
    code = "".join(code)
    
    mainBody = (
                    "\"(patch1|patch2|patch3)\"                     \n"     
                    "{                                              \n"     
                    "\t// Neumann boundary                          \n"     
                    "\ttype            codedMixed;                  \n"      
                    "\trefValue        uniform %s;                  \n"
                    "\trefGradient     uniform %s;                  \n"
                    "\tvalueFraction   uniform 0;                   \n"
                    "                                               \n" 
                    "\tname        %s;                              \n"  
                    "                                               \n" 
                    "\tcode                                         \n" 
                    "\t#{                                           \n" 
                    "                                               \n"
                    "\t\t// Gets current patch                      \n"
                    "\t\tconst fvPatch& boundaryPatch = patch();    \n"
                    "                                               \n"
                    "\t\t// Gets the patch face centres values      \n"
                    "\t\tconst vectorField& Cf = boundaryPatch.Cf(); \n"
                    "                                               \n"
                    "\t\tconst vectorField nf = patch().nf();       \n"
                    "                                               \n"
                    "\t\t// MMS                                     \n"
                    "                                               \n"
                    "%s                                             \n"
                    "\t\t// Loops over the patch                    \n"
                    "\t\tforAll(this->patch(), faceI)               \n"
                    "\t\t{                                          \n"
                    "%s                                             \n"
                    "%s                                             \n"
                    "%s                                             \n"
                    "\t\t}                                          \n\n"
                    "\t#};                                          \n"
                    "}"
                )
        
    print(mainBody
          %(
              value,
              value,
              name,
              transient.format("\t\t"),
              dimensions.format("\t\t\t", "Cf","[faceI]"),
              code.format("\t\t\t"),
              solution
              )
         )
    



def generateFunctionObject(MMS, variableName):
    isScalar = isinstance(MMS, sym.Expr)
    isVector = ( isinstance(MMS, sym.matrices.dense.MutableDenseMatrix) and MMS.shape == (3, 1) )

   
    transient = getTransient(MMS, True)
    dimensions = getDimensions(MMS)
    
    typeOfField = getTypeOfField(isScalar, isVector)
    typeOfField = "vol" + typeOfField[0].upper() + typeOfField[1:]
    
    mag = "mag"

    MMSFieldName = "MMS_diff_" + variableName
   
    if(isScalar):
        errorNorms = (
                          "Info << \"L1 norm is: \"    << gSum( {0}*V )/gSum(V)            << endl;     \n\n"  
                          "{1}Info << \"L2 norm is: \"    << sqrt( gSum({0}*{0}*V)/gSum(V) )    << endl;     \n\n"               
                          "{1}Info << \"Linf norm is: \"  << gMax( {0} )                 << endl;     \n\n"
                      ).format(MMSFieldName,"\t\t\t")

    if(isVector):
        errorNorms = (
                          "Info << \"For the 1st component of the vector\"                               << endl;     \n\n"
                          "{1}Info << \"L1 norm is: \"    << gSum( {0}.component(0)*V )/gSum(V)     << endl;     \n\n"  
                          "{1}Info << \"L2 norm is: \"    << sqrt( gSum( {0}.component(0)*{0}.component(0)*V)/gSum(V) ) << endl;     \n\n"               
                          "{1}Info << \"Linf norm is: \"  << gMax( {0}.component(0).ref() )             << endl;     \n\n\n"
                          
                          "{1}Info << \"For the 2nd component of the vector\"                            << endl;     \n\n"
                          "{1}Info << \"L1 norm is: \"    << gSum( {0}.component(1)*V )/gSum(V)     << endl;     \n\n"  
                          "{1}Info << \"L2 norm is: \"    << sqrt( gSum( {0}.component(1)*{0}.component(1)*V)/gSum(V) ) << endl;     \n\n"               
                          "{1}Info << \"Linf norm is: \"  << gMax( {0}.component(1).ref() )             << endl;     \n\n\n"
                          
                          "{1}Info << \"For the 3rd component of the vector\"                            << endl;     \n\n"
                          "{1}Info << \"L1 norm is: \"    << gSum( {0}.component(2)*V )/gSum(V)     << endl;     \n\n"  
                          "{1}Info << \"L2 norm is: \"    << sqrt( gSum( {0}.component(2)*{0}.component(2)*V)/gSum(V) ) << endl;     \n\n"               
                          "{1}Info << \"Linf norm is: \"  << gMax( {0}.component(2).ref())             << endl;     \n\n"
                      ).format(MMSFieldName,"\t\t\t")
        
    
    vectorSolution = "{0}"
    field = (
                "{0}(                                \n"
                "{0}\tIOobject                       \n"
                "{0}\t(                              \n"
                "{0}\t\t\"{1}\",                \n"
                "{0}\t\tmesh().time().timeName(),    \n"
                "{0}\t\tmesh(),                      \n"
                "{0}\t\tIOobject::NO_READ,           \n"
                "{0}\t\tIOobject::AUTO_WRITE         \n"
                "{0}\t),                             \n"
                "{0}\tmesh(),                         \n"
                "{0}\tdimensionedScalar (\"{1}_\", dimless, 0.0) \n"
                "{0});                               \n\n"
            )   
    
    if(isVector):
        vectorSolution = "{0}const vector solution (solution_1, solution_2, solution_3); \n"
        field = (
            "{0}(                                \n"
            "{0}\tIOobject                       \n"
            "{0}\t(                              \n"
            "{0}\t\t\"{1}\",                     \n"
            "{0}\t\tmesh().time().timeName(),    \n"
            "{0}\t\tmesh(),                      \n"
            "{0}\t\tIOobject::NO_READ,           \n"
            "{0}\t\tIOobject::AUTO_WRITE         \n"
            "{0}\t),                             \n"
            "{0}\tmesh(),                         \n"
            "{0}\tdimensionedVector (\"{1}_\", dimless, vector(0, 0, 0) ) \n"
            "{0});                               \n\n"
        )
        mag="cmptMag"
        
    mainBody = (
                    "functions                              \n"
                    "{{                                     \n"
                	"\terrorNorm_{2}                           \n"
                	"\t{{                                  \n"
                	"\t\ttype coded;                       \n"
                	"\t\tlibs (utilityFunctionObjects);    \n"
                	"\t\twriteControl writeTime;           \n"
                    "                                      \n"
                	"\t\tname analyticalSolution_{2};          \n"
                	"                                      \n"    
                	"\t\tcodeWrite                         \n"
                	"\t\t#{{                               \n"
                    "\t\t\tconst {0}& {1} = mesh().lookupObject<{0}>(\"{2}\"); \n\n"
                    "\t\t\tconst volVectorField& C = mesh().C(); \n\n"
                    "\t\t\tconst surfaceVectorField& Cf = mesh().Cf(); \n\n"
                    "\t\t\tconst scalarField& V = mesh().V();   \n\n"
                    "\t\t\t{0} {11}                     \n"
                    "{3}"

                    "{4}                                    \n\n"
                    "\t\t\tforAll({11}, cellI)          \n"
                    "\t\t\t{{                               \n"
                    "{5}                                    \n"
                    "{6}                                    \n"
                    "{7}                                    \n"
                    "\t\t\t\t{11}[cellI] = {10}(solution - {1}[cellI]); \n"
                    "\t\t\t}}                               \n\n"
                    
                    "\t\t\tforAll({11}.boundaryField(), patchI)      \n"
                    "\t\t\t{{                               \n"
                    "\t\t\t\tforAll({11}.boundaryField()[patchI], faceI)      \n"
                    "\t\t\t\t{{                             \n"
                    "{8}                                    \n"
                    "\t{6}                                  \n"
                    "\t{7}                                  \n"
                    "\t\t\t\t\t{11}.boundaryFieldRef()[patchI][faceI] = {10}(solution - {1}.boundaryField()[patchI][faceI]); \n"
                    "\t\t\t\t}}                             \n"
                    "\t\t\t}}                               \n\n"
                    "\t\t\t{9}                              \n"
                    "\t\t\t{11}.write();                \n"
                	"\t\t#}};                               \n"
                	"\t}}                                   \n"
                    "}}"    
                 ).format(
                             typeOfField,                       
                             variableName + "_find",
                             variableName,
                             field.format("\t\t\t", MMSFieldName),
                             transient.format("\t\t\t"),
                             dimensions.format("\t\t\t\t", "C", "[cellI]"),
                             generateCcode(MMS, "solution").format("\t\t\t\t"),
                             vectorSolution.format("\t\t\t\t"),
                             dimensions.format("\t\t\t\t\t", "Cf.boundaryField()", "[patchI][faceI]"),
                             errorNorms,
                             mag,
                             MMSFieldName
                             )
                    
    print(mainBody)
    
    

    



def generateFvOptions( sourceTerm,
                       titleForSource,
                       variableInWhichToApllySourceTerm,                
                     ):
    
    if(isinstance(sourceTerm, sym.matrices.immutable.ImmutableDenseMatrix)):
       sourceTerm=sourceTerm.as_mutable()
    
    isScalar = isinstance(sourceTerm, sym.Expr)
    isVector = ( isinstance(sourceTerm, sym.matrices.dense.MutableDenseMatrix) and sourceTerm.shape == (3, 1) )

    vectorOrScalarSourceCode = getTypeOfField(isScalar, isVector).strip("Field")
    
    transient= getTransient(sourceTerm, fvOptions=True)
    
    dimensions = getDimensions(sourceTerm)
    
    solution = ""
    
    if (isScalar):
        solution = ("\t\t\t{0}[cellI] -= V[cellI]*({1});").format(variableInWhichToApllySourceTerm+"Source", "solution_"+variableInWhichToApllySourceTerm)
        
    elif (isVector):
        solution = ( "\n\t\t\tconst vector solution ("+ ", ".join(standardSyntax("solution_" + variableInWhichToApllySourceTerm)) + "); \n"
    
                "\t\t\t{0}[cellI] -= V[cellI]*solution; \n" 
              ).format(variableInWhichToApllySourceTerm+"Source")
    else:
        raise ValueError('variable vectorOrScalarSourceCode is not defined correctly. \n Should be either scalar or vector')
    
    
    formatedCode = generateCcode(sourceTerm, "solution_" + variableInWhichToApllySourceTerm)
    
    Body = (
            "/*--------------------------------*- C++ -*----------------------------------*\ \n"
            "| =========                 |                                                |  \n"
            "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n"
            "|  \\    /   O peration     | Version:  v2012                                 | \n"
            "|   \\  /    A nd           | Website:  www.openfoam.com                      | \n"
            "|    \\/     M anipulation  |                                                 | \n"
            "\*---------------------------------------------------------------------------*/ \n"
            "FoamFile                                                                        \n"
            "{                                                                               \n"
            "    version     2.0;                                                            \n"
            "    format      ascii;                                                          \n"
            "    class       dictionary;                                                     \n"
            "    object      fvOptions;                                                      \n"
            "}                                                                               \n"
            "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n"
            "\n%s                                                                            \n"
            "{                                                                               \n"                                                                                
            "\ttype            %s;                                                           \n"
            "\tselectionMode   all;                                                          \n"
            "\tfields          (%s);                                                         \n"
            "                                                                                \n"                                                                              
            "\t// Name of the coded source                                                   \n"
            "\tname            %s;                                                           \n"
            "                                                                                \n"                                                                                
            "\tcodeInclude                                                                   \n"
            "\t#{                                                                            \n"
            "\t\t// Info: Include necessary libraries for calculation                        \n"
            "\t#};                                                                           \n"
            "                                                                                \n"                                                                                
            "\tcodeCorrect                                                                   \n"
            "\t#{                                                                            \n"
            "\t\t// Info: Apply corrections after the equation has been solved               \n"
            "\t#};                                                                           \n"
            "                                                                                \n"                                                                             
            "\tcodeConstrain                                                                 \n"
            "\t#{                                                                            \n"
            "\t\t// Info: Constrain values before the equation is solved                     \n"
            "\t#};                                                                           \n"                                                                                                                                                         
            )

    source = (
            "\tcodeAddSup                                                                   \n"
            "\t#{{                                                                          \n"                                                                        
            "\t\t// Gets the cell volumes of the mesh                                       \n"
            "\t\tconst scalarField& V = mesh_.V();                                          \n"
            "                                                                               \n"                                                                                  
            "\t\t// Gets the vector containing cell center position of the mesh             \n"
            "\t\tconst volVectorField& C = mesh().C();                                      \n"
            "                                                                               \n"                                                                                  
            "\t\t// Gets the equation source term                                           \n"
            "\t\t{0}& {1}= eqn.source();                                                    \n"
            "{2}                                                                             \n"
            "\t\t// Loops over each cell in the domain                                      \n"
            "\t\tforAll({1}, cellI)                                                         \n"
            "\t\t{{                                                                         \n"
            "{3}                                                                            \n"
            "{4}                                                                            \n"
            "{5}                                                                            \n"
            "\t\t}};                                                                         \n"                                                                                
            "\t#}};                                                                         \n"
            "}}                                                                             \n"
            ).format(
                        vectorOrScalarSourceCode + "Field",
                        variableInWhichToApllySourceTerm + "Source",
                        transient.format("\t\t"),
                        dimensions.format("\t\t\t", "C", "[cellI]"),
                        formatedCode.format("\t\t\t"),
                        solution
                    )


    print(
          Body 
          %(
              titleForSource,
              vectorOrScalarSourceCode + "CodedSource",
              variableInWhichToApllySourceTerm,
              titleForSource+"_2"
              )
          )
     
    print(source)
    
    
    

    
    

