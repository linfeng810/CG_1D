# CG_1D

This is an exercise implementation of continuous Galerkin solver for 1D unsteady diffusion problem. 

Modal problem:
\[
\partial p / \partial t - \partial^2 p / \partial x^2 = f(x,t) in \Omega \times [0,t_{end}) \\
p(0) = 0 in {0} \times [0,t_{end}) \\
p(1) = 0 in {0} \times [0,t_{end}) \\
p(x) = 0 in \Omega \times {0}
]\

Element: linear CG element

Quadrature: 1st order 2-point

Time discretisation: Backward Euler (fully implicit)
