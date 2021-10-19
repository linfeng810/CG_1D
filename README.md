# CG_1D

This is an exercise implementation of continuous Galerkin solver for 1D unsteady diffusion problem. 

Modal problem:
- <img src="https://latex.codecogs.com/gif.latex?\partial p / \partial t - \partial^2 p / \partial x^2 = f(x,t) \quad \text{in} \Omega \times [0,t_{end})" /> 
- <img src="https://latex.codecogs.com/gif.latex?p = 0 \quad \text{on } \partial\Omega \times [0,t_{end})" /> 
- <img src="https://latex.codecogs.com/gif.latex?p(x) = 0 in \Omega \times {0}" /> 

Element: linear CG element

Quadrature: 1st order 2-point

Time discretisation: Backward Euler (fully implicit)
