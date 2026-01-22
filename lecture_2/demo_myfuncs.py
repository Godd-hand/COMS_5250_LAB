# Python Lab Assignment

# 1. Factorial
"""
Aim: to calculate n!
Attempt: use a simple for loop
"""
def factorial(n):
    s = 1
    for k in range(1, n+1):
        s = s*(k)
    return s


# 2. Exponential Function
"""
Aim: to calculate the exponential function e^x.

Class Notes:

We will store the value of e ≈ 2.7182818284590451.
From x we find the nearest integer, let’s call it x_0:
x_0 = int(round(x))
We then compute the Taylor series of e^x about x = x_0.
"""
def exponential(x):
    e = 2.7182818284590451
    x0 = int(round(x))
    term = 1.0  # First term of the Taylor series
    sum_exp = term  # Initialize sum of series with the first term
    n = 1  # Counter for factorial in denominator

    while True:
        term = term * (x - x0) / n  # Compute next term in series
        sum_exp += term  # Add the new term to the sum
        if abs(term) < 1e-10:  # Convergence check
            break
        n += 1

    return sum_exp * (e ** x0)  # Scale by e^x0



# 3. Natural Logarithm Function
"""
Aim: to calculate the natural logarithm ln(X)
s = ln(x) 
e^s = x
f(s) = e^s - x
"""
def natural_log(x):
    if x <= 0:
        raise ValueError("Input must be a positive number.")
    
    # Initial guess for s
    s = x if x < 2 else 1.0
    tolerance = 1e-10
    max_iterations = 1000
    iteration = 0

    while iteration < max_iterations:
        f_s = exponential(s) - x
        f_prime_s = exponential(s)  # Derivative of e^s is e^s

        # Newton's method update
        s_new = s - f_s / f_prime_s

        # Check for convergence
        if abs(s_new - s) < tolerance:
            return s_new
        
        s = s_new
        iteration += 1

    raise RuntimeError("Failed to converge to a solution.")

