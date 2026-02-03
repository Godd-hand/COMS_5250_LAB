# Develop a function GaussElimination(A,b) for solving linear systems Ax = b
# by Gaussian Elimination, with A and b as inputs and x as output.

import numpy as np
def GaussElimination(A, b):
    """
    This function solves the linear system Ax = b using Gaussian Elimination and back substitution.

    Parameters:
    A is an ndarray of shape (n,n) representing the coefficient matrix a_11 to a_nn.
    b is an ndarray of shape (n,1) representing the right-hand side vector.

    Returns:
    x is an ndarray of shape (n,1) representing the solution vector.
    """
    n = A.shape[0] # to get the number of equations
    A = A.astype(float)  # ensure A is of type float for division operations
    b = b.astype(float)  # ensure b is of type float for division operations

    # Build augmented matrix [A | b]
    M = np.zeros((n, n + 1), dtype=float)
    for i in range(n):
        for j in range(n):
            M[i, j] = A[i, j]
        M[i, n] = b[i, 0]
    print("Initial augmented matrix:")
    print(M)

    tol = 1e-12
    for i in range(n):
        max_row = i + np.argmax(np.abs(M[i:, i])) # partial pivoting with i added to index correctly
        if max_row != i:
            M[[i, max_row], :] = M[[max_row, i], :] # swap rows if needed 

        pivot = M[i, i]
        if abs(pivot) < tol:
            raise ValueError("Matrix is singular or nearly singular; cannot proceed with Gaussian elimination")

        for j in range(i + 1, n):
            factor = M[j, i] / pivot
            for k in range(i, n + 1):
                M[j, k] = M[j, k] - factor * M[i, k]

    print("Augmented matrix after forward elimination:")
    print(M)
    x = np.zeros(n) # Initialize the solution vector
    for i in range(n - 1, -1, -1): # Start from the last row and move upwards
        if abs(M[i, i]) < tol:
            raise ValueError("Zero pivot encountered during back substitution; system is singular or ill-conditioned")
        s = 0.0
        for j in range(i + 1, n):
            s += M[i, j] * x[j]
        x[i] = (M[i, n] - s) / M[i, i]
        #print(x)
    
    return x.reshape((n, 1))


xs = [-0.1, -0.02, 0.02, 0.1]
A = np.zeros((4, 4), dtype=float) # Initialize coefficient matrix
b = np.zeros((4, 1), dtype=float) # Initialize right-hand side vector
for i in range(4):
    x = xs[i]
    A[i, 0] = x**3
    A[i, 1] = x**2
    A[i, 2] = x
    A[i, 3] = 1.0
    b[i, 0] = np.cos(x)
print("Coefficient matrix A:")
print(A)
print("Right-hand side vector b:")
print(b)

coefficients = GaussElimination(A, b)

a, b_coef, c, d = coefficients.flatten()
print("Polynomial coefficients:")
print(f"a = {a}, b = {b_coef}, c = {c}, d = {d}")

def p(x):
    return a * x**3 + b_coef * x**2 + c * x + d

max_error = 0.0
for x in xs:
    err = abs(np.cos(x) - p(x))
    if err > max_error:
        max_error = err
print(f"Maximum absolute error: {max_error}")
