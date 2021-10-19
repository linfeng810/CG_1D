# CG_1D

This is an exercise implementation of continuous Galerkin solver for 1D unsteady diffusion problem. 

Modal problem:

<a href="https://www.codecogs.com/eqnedit.php?latex=\partial&space;p&space;/&space;\partial&space;t&space;-&space;\partial^2&space;p&space;/&space;\partial&space;x^2&space;=&space;f(x,t)&space;\quad&space;\text{in&space;}&space;\Omega&space;\times&space;[0,t_{end})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\partial&space;p&space;/&space;\partial&space;t&space;-&space;\partial^2&space;p&space;/&space;\partial&space;x^2&space;=&space;f(x,t)&space;\quad&space;\text{in&space;}&space;\Omega&space;\times&space;[0,t_{end})" title="\partial p / \partial t - \partial^2 p / \partial x^2 = f(x,t) \quad \text{in } \Omega \times [0,t_{end})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=p&space;=&space;0&space;\quad&space;\text{on&space;}&space;\partial\Omega&space;\times&space;[0,t_{end})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p&space;=&space;0&space;\quad&space;\text{on&space;}&space;\partial\Omega&space;\times&space;[0,t_{end})" title="p = 0 \quad \text{on } \partial\Omega \times [0,t_{end})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=p(x)&space;=&space;0&space;\quad\text{in&space;}&space;\Omega&space;\times&space;\{0\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p(x)&space;=&space;0&space;\quad\text{in&space;}&space;\Omega&space;\times&space;\{0\}" title="p(x) = 0 \quad\text{in } \Omega \times \{0\}" /></a>

Element: linear CG element

Quadrature: 1st order 2-point

Time discretisation: Backward Euler (fully implicit)
