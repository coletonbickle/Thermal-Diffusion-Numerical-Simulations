# Thermal-Diffusion-Numerical-Simulations


Description: One Dimensional Thermal Diffusion Problem

Author: Coleton C. Bickle

Methods:
    Pseudo-Transient Steady State Solver,
    Second Order Accurate Spatial Discretization,
    Gauss-Seidel Matrix-Wise Linear Solver,
    and an Under Relaxation Factor
    
One Dimensional Boundary Conditions:
    Isothermal (Left Side)
    Convective Heat Transfer (Right Side)
    Source Term (Heat Generation Within Rod)

This One Dimensional Thermal Diffusion Problem uses a Pseudo-Transient Steady-State solver while implementing an Isothermal, Heat convection, and Source Term boundary conditions.
These implementations display a good non-linearity in temperature change as well as a temperature drop from the right side of the rod due to the heat convection boundary condition.


![1d_solution2](https://github.com/coletonbickle/Thermal-Diffusion-Numerical-Simulations/assets/91445808/71132ea6-adbe-44c2-b524-9bb1e1d0f4d4)

Description: Two Dimensional Thermal Diffusion Problem

Author: Coleton C. Bickle

Methods:
    Pseudo-Transient Steady State Solver,
    Second Order Accurate Spatial Discretization,
    Gauss-Seidel Matrix-Wise Linear Solver,
    Under Relaxation Factor,
    
Two Dimensional Boundary Conditions:
    Isothermal Wall,
    Convective Heat Transfer,
    Adiabatic Wall,
    Heat Sink,
    Source Term,

This Two Dimensial Thermal Diffusion Problem is an expansion on the 1D problem around the same Heating Equation. This problem simulates a copper plate with various boundary conditions applied to each edge. At the same time, there is heat generation source term applied evenly throughout the plate.

![2d_solution](https://github.com/coletonbickle/Thermal-Diffusion-Numerical-Simulations/assets/91445808/24729289-892f-448d-b6a5-d4f38c2e2d0c)


Future Work:
    After learning a lot regarding Implicit Steady State Solvers, I plan on
    attempting to code the 1D Sod Shock Tube problem to dive into the realm
    of compressible flows.


