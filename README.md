_pyMMSFOAM_ is a python module for generating manufactured solutions for verification of solvers implemented in OpenFOAM&reg;.

The utility was developed based on OpenFOAM&reg; V2106 and resorts to the _coded_ framework to implement _boundaryConditions_, _functionObjects_ and _fvOptions_.

_pyMMSFOAM_ is based on [Sympy](https://www.sympy.org/en/index.html) and is composed of different submodules, the most important being:
* Symbols:
    * Represents the symbols that the user is allowed to use for symbolic differentiation. Currently it only allows the usage of symbols representing the coordinates of Cartesian reference frame, i.e., $x$ , $y$, and $z$.
* Differentiation operators:
    * Has functions to compute the symbolic differentiation operators tipically used in the finite volume method:
        * Gradient, $\nabla$,
        * Divergence, $\nabla \cdot$,
        * Laplacian,  $\nabla \left(\Gamma \nabla \cdot \right)$
        * If required other operators can be easily implemented in _pyMMSFOAM_.
* Printing utilities:
    * Is responsible for printing the source terms resulting from the method of manufactured solutions, suitable _boundary conditions_ (fixed value or fixed gradient) as well as _functionObjects_ for computing several common error norms ($L^1$, $L^2$, $L^\infty$).
    * The output should be copied into the respective case study locations.

This repository is associated to a paper published in the OpenFOAM Journal [1], which describes the case studies named _laplacianFoam.py_ and _simpleFoam.py_, and should provide suitable tutorials for using the python utility.

The authors welcome suggestions for improvements as well as manufactured solutions for creating new case study suits.

## References
[1] - B. Ramoa et al. "A semi-automatic approach based on the method of manufactured solutions to assess the convergence order in OpenFOAM" , OpenFOAM Journal, <add pages when available>, 2022, <doi when available>.
