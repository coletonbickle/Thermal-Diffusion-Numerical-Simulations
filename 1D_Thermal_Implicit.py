import numpy as np
import matplotlib.pyplot as plt

'''
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
'''


def main():
    
    # Spatial Discretization Properties
    pts = 50 # Total Grid Points
    vnn = 4  # Von Neumann Number
    L = 1    # Length Scale[m]
    
    # Pseudo-Time and Spatial Control
    nt = 100000 # Max allowable iterations
    itol = 1e-5 # Simulation residual tolerance
    gtol = 0.1  # G-S tolerance
    gmax = 50   # Max G-S iterations
    URF = 0.7   # Under Relaxation Factor
    
        
    # Material Properties (Copper)
    alpha = 1.17e-4 # Thermal Diffusivity [m^2/s]
    cp = 385        # Specific Heat Capacity [J/kg*K]
    rho = 8960      # Density [Kg/m^3]
    
    # Initial Conditions
    T_L = 273  # Initial Temp Left side of rod [K]
    T_R = 323  # Initial Temp right side of rod [K]
    T_IC = 323 # Initial Temp inside rod [K]
    
    # Boundary Conditions
    qdot = 10   # Internal heating (source term)
    h = 100     # Convective Coefficient of external fluid (Water)
    T_inf = 280 # Ambient Temperature of external fluid (Water)
    

    # Pre-Process
    k = alpha*rho*cp       # Thermal Conduction of Copper
    S = qdot*1000/rho/cp   # Heat Source Term
    dx = L/(pts-1)         # Distance between nodes
    dtau = vnn*dx*dx/alpha # Pseudo Time Step
    Stau = dtau*S          # Discretized Source Term
    Nu = h*dx/k            # Nusselt Number
    Nuinv = 1/Nu           # Set up 1/Nu for Heat Convection BC
    
    # Set up Initial Conditions
    bi = np.ones((pts,1))*T_IC
    bi[0] = T_L
    bi[pts-1:] = T_R
    b = bi.copy()
    
    # Generate A Matrix
    A = create_A(pts,vnn)

    # Set up Matrix for G-S Solver
    U = np.triu(A,1)
    L = np.tril(A)
    LI = np.linalg.inv(L)
    T = -np.dot(LI,U)
    
    # Solver
    i = 1
    ierr = 1
    while (nt > i) and (ierr > itol):
        # Gauss-Seidel Matrix-Wise Method
        x = GS_matrix(T,LI,b,gtol,gmax)
        
        # Apply Boundary Conditions
        x[1:-2] = x[1:-2] + Stau
        x[-1] = (T_inf + Nuinv*x[-2])/(Nuinv+1)
        
        # Calculate Residual
        ierr = np.max(np.sqrt((b-x)*(b-x)))
        
        # Apply Under Relaxation Factor
        b = x.copy()*URF + (1-URF)*b.copy()

        # Print Progress
        if i % 100 == 0:
            print(f"Iteration: %d" % i,"Residual:","{:e}".format(ierr))
        i = i + 1


    # Post-Process
    print("\n==============DONE==============")
    print(f"\nTotal Iterations: %d" % i,"\nFinal Residual:","{:e}".format(ierr),"\n")
    
    X = np.linspace(0,pts,pts)
    plt.figure()
    plt.plot(X,b)
    plt.plot(X,bi)
    plt.xlabel("Node Number")
    plt.ylabel("Temperature [K]")
    plt.title("1D Thermal Diffusion Problem")
    plt.xlim(0,pts)
    plt.legend(['Solution','Initial Condition'])
    plt.show()
    
def create_A(pts,vnn):
    A = np.zeros((pts,pts))

    for i in range(pts):
        for j in range(pts):
            if i == j:
                A[i,j] = (1+2*vnn)
            if j == (i+1):
                A[i,j] = (-1*vnn)
            if j == (i-1):
                A[i,j] = (-1*vnn)
    A[0,0] = 1
    A[0,1] = 0
    A[-1,-1] = 1
    A[-1,-2] = 0
    
    return A

def GS_matrix(T,LI,b,gtol,gmax):
    k = 1
    gerr = 1
    x = b.copy()
    C = np.dot(LI,x)
    while (np.max(gerr) > gtol) and (k < gmax):
        xp = np.dot(T,x) + C
        gerr = np.sqrt((x-xp)*(x-xp))
        x = xp
        k = k + 1
    
    return x

if __name__ == '__main__':
    main()