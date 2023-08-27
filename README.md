# Thermal-Diffusion-Numerical-Simulations


Description: One Dimensional Thermal Diffusion Problem

Author: Coleton C. Bickle

Methods:
    Pseudo-Transient Second Order Accurate Spatial Discretization
    Gauss-Seidel Matrix-Wise Linear Solver
    Under Relaxation Factor
    
Boundary Conditions
    Isothermal (Left Side)
    Convective Heat Transfer (Right Side)
    Source Term (Heat Generation Within Rod)
    
Future Work:
    I am intending to apply what I have learned from implementing an implicit
    numerical method to a one dimensional thermal diffusion problem eventually
    to a two dimensional thermal diffusion across a plate while implementing
    various different boundary conditions

This One Dimensional Thermal Diffusion Problem uses Pseudo-Transient solution marching while implementing an Isothermal, Heat convection, and Source Term boundary conditions
These implementations display a good non-linearity in temperature change as well as a temperature drop from the right side of the rod due to the heat convection boundary condition
