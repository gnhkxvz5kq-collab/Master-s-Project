import numpy as np
from scipy.sparse import diags
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import spsolve

def cell_locations(no_cells, cell_length):
    """
    Return total fibre length and cell-centre positions.
    """
    length = no_cells * cell_length
    positions = (np.arange(no_cells) + 0.5) * cell_length

    return length, positions

def spatial_derivative_1D(N, L):
    """
    Construct the 1D second spatial derivative matrix
    using a cell-centred finite-volume discretisation.

    Parameters
    ----------
    N : int
        Number of spatial control volumes.
    L : float
        Total length of the domain (m).
 
    Returns
    -------
    D2 : scipy.sparse.csr_matrix
        Sparse second derivative matrix (1/m^2).
    x : numpy.ndarray
        Positions of the cell centres (m).
    dx : float
        Spatial grid spacing (m).
    """

    if not isinstance(N, (int, np.integer)) or N < 1:
        raise ValueError("N must be a positive integer.")

    if L <= 0:
        raise ValueError("L must be positive.")

    dx = L / N

    # Cell-centre positions
    x = (np.arange(N) + 0.5) * dx

    # Main diagonal
    main = -2.0 * np.ones(N)

    # Zero-flux boundary conditions
    main[0] = -1.0
    main[-1] = -1.0

    # Adjacent diagonals
    upper = np.ones(N - 1)
    lower = np.ones(N - 1)

    # Construct sparse matrix
    D2 = diags(
        [lower, main, upper],
        offsets=[-1, 0, 1],
        shape=(N, N),
        format="csr"
    ) / dx**2

    return D2, x, dx

'''below function defines the conductivity matrices if working only in 1D - these are isotropic so sigma_ie are both constants
'''
def conductivity_matrices(D2, sigma_i, sigma_e):
    """
    Construct intracellular and extracellular
    conductivity matrices.

    Parameters
    ----------
    D2 : sparse matrix
        Second spatial derivative matrix.
    sigma_i : float
        Intracellular conductivity (S/m).
    sigma_e : float
        Extracellular conductivity (S/m).

    Returns
    -------
    Ki, Ke : sparse matrices
        Conductivity stiffness matrices.
    """

    if sigma_i <= 0 or sigma_e <= 0:
        raise ValueError("Conductivities must be positive.")

    Ki = -sigma_i * D2
    Ke = -sigma_e * D2

    return Ki, Ke

'''given matrix equation for derivatives of phi_e and V_m we can find the inverse matrices to solve this for phi_e / V_m
need to define function to invert the matrix
this function uses a built in function from scipy to solve the matrix equation
'''
def solve_extracellular_potential(Ki, Ke, Vm):
    """
    Solve the global extracellular potential equation:

        (Ki + Ke) @ phi_e = -Ki @ Vm

    with a zero-mean extracellular potential constraint.

    Parameters
    ----------
    Ki : sparse matrix
        Intracellular conductivity stiffness matrix.
    Ke : sparse matrix
        Extracellular conductivity stiffness matrix.
    Vm : numpy array
        Membrane potential at each spatial point.

    Returns
    -------
    phi_e : numpy array
        Extracellular potential at each spatial point.
    """

    N = len(Vm)

    # Construct global extracellular matrix
    A = Ki + Ke

    # Construct right-hand side
    b = -(Ki @ Vm)

    # Add zero-mean gauge constraint
    ones = csr_matrix(np.ones((N, 1)))

    A_aug = bmat([
        [A, ones],
        [ones.T, None]
    ], format="csr")

    b_aug = np.append(b, 0.0)

    # Solve the augmented sparse linear system
    solution = spsolve(A_aug, b_aug)

    # Extract extracellular potential
    phi_e = solution[:N]

    return phi_e

'''now we have ways to solve for phi_e define an equation to solve for phi_i'''

def solve_phi_i(phi_e, V_m):
    phi_i = V_m + phi_e
    return phi_i

'''calculate from the differential equation the rate of change of voltage w.r.t time'''

def calculate_dVm(Ki, Vm, phi_e, beta, Cm,
                  I_ion, I_stim):
    """
    Calculate the rate of change of membrane voltage.

    Assumes:
        - Constant beta and Cm
        - Consistent spatial units
        - Outward ionic current is positive
        - Inward stimulus current is positive
    """

    phi_i = Vm + phi_e

    dVm_dt = (
        -(Ki @ phi_i) / (beta * Cm)
        - I_ion / Cm
        + I_stim / Cm
    )

    return dVm_dt

'''this dv/dt can then be used in euler method integration to find the membrane potential at the next time interval
the new values for phi_i, phi_e, and V_m can then be used to repeat this process
Now we have the new voltage, this will be used in the Purkinje model as an input to calculate the current flows
and the new state variables 
This new voltage is also used in the bidomain functions to calculate the new intra and extracellular potentials and 
new dV/dt to be used for the next time step
'''
'''calculate the new membrane voltage'''
def membrane_potential(Vm, dVm_dt, dt):
    Vm = + dt * dVm_dt
    return Vm


'''define a bidomain function that combines all above elements'''

def bidomain(Vm, no_cells, cell_length, sigma_i, sigma_e, I_ion, I_stim, beta, Cm, dt):
    length, positions = cell_locations(no_cells, cell_length)
    D2, x, dx = spatial_derivative_1D(no_cells, length)
    Ki, Ke = conductivity_matrices(D2, sigma_i, sigma_e)
    phi_e = solve_extracellular_potential(Ki, Ke, Vm)
    phi_i = Vm + phi_e
    dVm_dt = calculate_dVm(Ki, Vm, phi_e, beta, Cm, I_ion, I_stim)
    Vm_new = membrane_potential(Vm, dVm_dt, dt)
    return Vm_new

