# Thermal-Diffusion-Numerical-Simulations


Description: One Dimensional Thermal Diffusion Problem

Author: Coleton C. Bickle

Methods:
    Pseudo-Transient Steady State Solver,
    Second Order Accurate Spatial Discretization,
    Gauss-Seidel Matrix-Wise Linear Solver,
    and an Under Relaxation Factor
    
Boundary Conditions:
    Isothermal (Left Side)
    Convective Heat Transfer (Right Side)
    Source Term (Heat Generation Within Rod)
    
Future Work:
    I am intending to apply what I have learned from implementing an implicit
    numerical method to a one dimensional thermal diffusion problem eventually
    to a two dimensional thermal diffusion problem across a heated plate while implementing
    various different boundary conditions

This One Dimensional Thermal Diffusion Problem uses a Pseudo-Transient Steady-State solver while implementing an Isothermal, Heat convection, and Source Term boundary conditions.
These implementations display a good non-linearity in temperature change as well as a temperature drop from the right side of the rod due to the heat convection boundary condition.


![1d_solution2](https://github.com/coletonbickle/Thermal-Diffusion-Numerical-Simulations/assets/91445808/71132ea6-adbe-44c2-b524-9bb1e1d0f4d4)



