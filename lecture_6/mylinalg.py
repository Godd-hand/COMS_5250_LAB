# Module: mylinalg.py

"""
This module contains the GaussElimination(A,b) function produced in the last assignment, together
with partial pivoting to improve numerical stability.

In addition, a new function LeastSquareApprox(x,f,n) is provided to find the least square approximation 
of data {x, f} by a polynomial of degree <= n. The output is the coefficient of the polynomial.

Application of this module:

1. Find the polynomial p(x) of degree <= 5 that approximates the function f(x) = cos(x) at the nodes linspace(-pi, pi, 51)
    in the least square sense.
2. Include the "main program" in this module as well.
3. Plot the function f(x) and the polynomial p(x) on the interval [-pi, pi].
4. Run the module as a script to produce the plot.

"""
import numpy as np
import matplotlib.pyplot as plt

# Gaussian Elimination function with partial pivoting
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
    # print("Initial augmented matrix:")
    # print(M)
    
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

    # print("Augmented matrix after forward elimination:")
    # print(M)
    x = np.zeros(n) # Initialize the solution vector
    for i in range(n - 1, -1, -1): # Start from the last row and move upwards
        if abs(M[i, i]) < tol:
            raise ValueError("Zero pivot encountered during back substitution; system is singular or ill-conditioned")
        s = 0.0
        for j in range(i + 1, n):
            s += M[i, j] * x[j]
        x[i] = (M[i, n] - s) / M[i, i]
        # print(x)
    
    return x.reshape((n, 1))

def LeastSquareApprox(x, f, n):
    """
    This function computes the least square approximation of data {x, f} by a polynomial of degree <= n.

    Parameters:
    x is an ndarray of shape (m,) representing the data points.
    f is an ndarray of shape (m,) representing the function values at the data points.
    n is an integer representing the degree of the approximating polynomial.

    Returns:
    a is an ndarray of shape (n+1, 1) representing the coefficients of the approximating polynomial.
    """
    m = len(x)
    A = np.zeros((n + 1, n + 1), dtype=float)
    b = np.zeros((n + 1, 1), dtype=float)

    for i in range(n + 1):
        for j in range(n + 1):
            A[i, j] = np.sum(x ** (i + j))
        b[i, 0] = np.sum(f * (x ** i))

    # print("Normal equation coefficient matrix A:")
    # print(A)
    # print("Normal equation right-hand side vector b:")
    # print(b)

    a = GaussElimination(A, b)
    return a

# Evaluating the polynomial at a given x
def EvaluatePolynomial(a, x):
    """
    This function evaluates the polynomial with coefficients a at the point x.

    Parameters:
    a is an ndarray of shape (n+1, 1) representing the coefficients of the polynomial.
    x is a float representing the point at which to evaluate the polynomial.

    Returns:
    p is a float representing the value of the polynomial at x.
    """
    n = a.shape[0] - 1
    p = 0.0
    for i in range(n + 1):
        p += a[i, 0] * (x ** i)
    return p


# Application: Find the least square polynomial approximation of f(x) = cos(x) on [-pi, pi]
if __name__ == "__main__":
    x = np.linspace(-np.pi, np.pi, 51)
    f = np.cos(x)
    n = 5
    a = LeastSquareApprox(x, f, n)
    print("Coefficients of the least square polynomial approximation:")
    print(a)

    # Plotting
    x_plot = np.linspace(-np.pi, np.pi, 100)
    f_plot = np.cos(x_plot)
    a_flat = a.flatten()
    a_flip = np.flip(a_flat)
    # print("Flipped coefficients for np.polyval:")
    # print(a_flip)
    p_plot = EvaluatePolynomial(a, x_plot)  # evaluate polynomial using coefficients
    # p_plot = np.polyval(a_flip, x_plot)  # evaluate polynomial using coefficients

    plt.figure(figsize=(10, 5))
    plt.plot(x_plot, f_plot, label='f(x) = cos(x)', linewidth=2)
    plt.plot(x_plot, p_plot, label=f'Least square polynomial of degree {n}', linewidth=2)
    plt.scatter(x, f, color='red', label='Data points')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Least Square Polynomial Approximation')
    plt.legend()
    plt.grid(True)
    plt.show()
