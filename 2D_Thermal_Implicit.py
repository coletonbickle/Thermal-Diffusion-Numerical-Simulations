import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

'''
Description: Two Dimensional Thermal Diffusion Problem
Author: Coleton C. Bickle

Methods:
    Pseudo-Transient Steady State Solver
    Second Order Accurate Spatial Discretization
    Gauss-Seidel Matrix-Wise Linear Solver
    Under Relaxation Factor
    
Boundary Conditions:
    Isothermal Wall
    Convective Heat Transfer
    Adiabatic Wall
    Heat Sink
    Source Term
    
Future Work:
    After learning a lot regarding Implicit Steady State Solvers, I plan on
    attempting to code the 1D Sod Shock Tube problem to dive into the realm
    of compressible flows.
'''


def main():
    
    # Spatial Discretization Properties
    pts = 55   # Total Grid Points
    vnn = 100  # Von Neumann Number
    L = 2      # Length Scale [m]
    
    # Pseudo-Time and Spatial Control
    nt = 100000  # Max allowable iterations
    itol = 1e-6  # Simulation residual tolerance
    gtol = 0.1   # G-S tolerance
    gmax = 30    # Max G-S iterations
    URF = 1      # Under Relaxation Factor
    
    # Material Properties (Copper)
    alpha = 1.17e-4 # Thermal Diffusivity [m^2/s]
    cp = 385        # Specific Heat Capacity [J/kg*K]
    rho = 8960      # Density [Kg/m^3]
    
    # Initial Conditions
    T_L = 300  # Initial Temp Left side of rod [K]
    T_R = 300  # Initial Temp Right side of rod [K]
    T_T = 300  # Initial Temp Top side of rod [K]
    T_B = 300  # Initial Temp Bottom side of rod [K]
    T_IC = 373 # Initial Temp inside rod [K]
    
    # Boundary Conditions
    qdot = 20     # Internal heating (source term) [W/m^2]
    q2dot = -8000 # Heat Sink [W/m^2]
    h = 100       # Convective Coefficient of external fluid (Water)
    T_inf = 273   # Ambient Temperature of external fluid (Water)
    

    # Pre-Process
    k1 = alpha*rho*cp      # Thermal Conduction of Copper
    S = qdot*1000/rho/cp   # Heat Source Term
    dx = L/(pts-1)         # Distance between nodes
    dtau = vnn*dx*dx/alpha # Pseudo Time Step
    Stau = dtau*S          # Discretized Source Term
    Nu = h*dx/k1           # Nusselt Number
    Nuinv = 1/Nu           # Set up 1/Nu for Heat Convection BC
    N2 = pts*pts           # A matrix dimensions

    # Generate A Matrix and b inital vector
    A, b, Flag_R, Flag_L, Flag_T, Flag_B = create_A(pts,N2,vnn,T_L,T_R,T_T,T_B,T_IC,Stau)
    Flag = np.concatenate((Flag_R, Flag_L, Flag_T, Flag_B),axis=None)
    
    # Set up Matrix for G-S Solver
    U = np.triu(A,1)
    L = np.tril(A)
    LI = np.linalg.inv(L)
    T = -np.dot(LI,U)
    
    # Implicit Solver
    i = 1
    ierr = 1
    bp = b.copy()
    x = b.copy()
    k = 0
    while (nt > i) and (ierr > itol):
        # Gauss-Seidel Matrix-Wise Method
        x = GS_matrix(T,LI,b,x,gtol,gmax)
        
        # Calculate Residual
        ierr = np.max(np.sqrt((bp-x)*(bp-x)))
        
        # Apply Under Relaxation Factor
        bp = x.copy()*URF + (1-URF)*b.copy()
        
        # Apply Boundary Conditions
        for k in range(N2):
            if k not in Flag:
                b[k] = bp[k] + Stau                      # Internal Source Term
            elif np.isin(k,Flag_R):
                b[k] = bp[k-1]                           # Adiabatic  Wall BC
            elif np.isin(k,Flag_L):
                b[k] = q2dot*dx/k1+bp[k+1]               # Isothermal Wall BC
            elif np.isin(k,Flag_B):
                b[k] = (T_inf+Nuinv*bp[k-pts])/(Nuinv+1) # Heat Convection BC
            elif np.isin(k,Flag_T):
                b[k] = bp[k]                             # Isothermal Wall BC
            
        # Print Progress
        if i % 100 == 0:
            print(f"Iteration: %d" % i,"Residual:","{:e}".format(ierr))
        i = i + 1


    # Post-Process
    print("\n==============DONE==============")
    print(f"\nTotal Iterations: %d" % i,"\nFinal Residual:","{:e}".format(ierr),"\n")
    
    
    # Reshape b solution into a matrix
    sol=np.reshape(x,(pts,pts))
    xmesh = np.linspace(0,pts,pts)
    X, Y = np.meshgrid(xmesh,xmesh)
    
    parulaColor = get_color()
    # Plot Contour
    fig = plt.figure(figsize=(13, 5), dpi=100)
    fig.add_subplot(1,2,1)
    CS = plt.contour(X, Y, sol,9,colors='k')#,9)#,cmap=cm.jet)
    plt.clabel(CS, inline=1,fmt=fmt,fontsize=10,colors='k')
    plt.contourf(X, Y, sol,9,cmap=parulaColor)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Temperature [K]',rotation=270,labelpad=15)
    plt.xlabel("Node Number")
    plt.ylabel("Node Number")
    plt.title("2D Filled Contour")
    
    ax = fig.add_subplot(1,2,2,projection='3d')
    surf = ax.plot_surface(X, Y, sol, cmap=parulaColor)#cm.viridis
    plt.title("3D Surface Representation")
    fig.suptitle("2D Thermal Diffusion Problem")
    plt.show()


def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} K"

'''
This Function generates the A matrix used in the linear Gauss-Seidel linear solver
It also flags the boundary nodes of the 2d mesh and correctly adjusts the corresponding
cell values in the A matrix to adjust for the boundary conditions.
'''
def create_A(pts,N2,vnn,T_L,T_R,T_T,T_B,T_IC,Stau):
    A = np.zeros((N2,N2))
    
    for i in range(N2):
        for j in range(N2):
            if i == j:
                A[i,j] = (1+4*vnn)
            if j == (i+1):
                A[i,j] = (-1*vnn)
            if j == (i-1):
                A[i,j] = (-1*vnn)
            if j == (i+pts):
                A[i,j] = (-1*vnn)
            if j == (i-pts):
                A[i,j] = (-1*vnn)
    # Set up Initial Conditions
    b = np.ones((N2,1))*T_IC + Stau
    
    # Create Boundary Flags
    Flag_R = np.zeros((1,pts-2))
    Flag_L = np.zeros((1,pts-2))
    Flag_T = np.zeros((1,pts))
    Flag_B = np.zeros((1,pts))

    # Flag boundary node values and apply initial conditions to initial b vector
    k = 0
    for i in range(pts):
        for j in range(pts):
            if (i == 0) and (j<(pts)):
                Flag_T[:,j] = k
                b[k] = T_B
            elif (k >= (N2-pts)):
                Flag_B[:,j] = k
                b[k] = T_T
            elif (np.mod(k,pts) == 0) and (k > (pts-1)) and (k < (N2-pts-1)):
                Flag_L[:,i-1] = k
                b[k] = T_L
            elif (np.mod(k,pts) == (pts-1)) and (k > (pts)) and (k <= (N2-pts)):
                Flag_R[:,i-1] = k
                b[k] = T_R
            k = k + 1
    
    # Use Flag values and adjust A matrix with boundary values
    Flag_concat = np.concatenate((Flag_R,Flag_L,Flag_T,Flag_B),axis=None)
    for i in range(N2):
        for j in range(N2):
            if (i == j) and np.isin(j,Flag_concat):
                A[i,:] = np.zeros((1,N2))
                A[i,j] = 1
    
    return A, b, Flag_R, Flag_L, Flag_T, Flag_B

'''
Gauss-Seidel Iterative Linear Solver
'''
def GS_matrix(T,LI,b,x,gtol,gmax):
    k = 1
    gerr = 1
    C = np.dot(LI,b)
    while (gtol < np.max(gerr)) and (k < gmax):
        xp = np.dot(T,x) + C
        gerr = np.sqrt((x-xp)*(x-xp))
        x = xp
        k = k + 1
    
    return x


'''
Custom contour color palate
'''
def get_color():
    cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905], 
    [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143], 
    [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952, 
    0.779247619], [0.1252714286, 0.3242428571, 0.8302714286], 
    [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238, 
    0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571], 
    [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571, 
    0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429], 
    [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667, 
    0.8467], [0.0779428571, 0.5039857143, 0.8383714286], 
    [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571, 
    0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429], 
    [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524, 
    0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048, 
    0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667], 
    [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381, 
    0.7607190476], [0.0383714286, 0.6742714286, 0.743552381], 
    [0.0589714286, 0.6837571429, 0.7253857143], 
    [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429], 
    [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429, 
    0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048], 
    [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619, 
    0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667], 
    [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524, 
    0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905], 
    [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476, 
    0.4493904762], [0.609852381, 0.7473142857, 0.4336857143], 
    [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333], 
    [0.7184095238, 0.7411333333, 0.3904761905], 
    [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667, 
    0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762], 
    [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217], 
    [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857, 
    0.2886428571], [0.9738952381, 0.7313952381, 0.266647619], 
    [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857, 
    0.2164142857], [0.9955333333, 0.7860571429, 0.196652381], 
    [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857], 
    [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309], 
    [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333, 
    0.0948380952], [0.9661, 0.9514428571, 0.0755333333], 
    [0.9763, 0.9831, 0.0538]]

    parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
    
    return parula_map



if __name__ == '__main__':
    main()