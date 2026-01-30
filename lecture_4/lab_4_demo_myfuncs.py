# Python Lab 4 Assignment

## Part 1: Update the functions, make them more robust to bad inputs and include checks

## constant e
e = 2.7182818284590451


## 1. Factorial
def factorial(n):
    """
    Compute n! for a non-negative integer n.

    Ensures `n` is an integer >= 0.
    """
    if not isinstance(n, int):
        raise TypeError("factorial() argument must be an integer") # updated to check type
    if n < 0:
        raise ValueError("factorial() not defined for negative values") # updated to check non-negativity
    
    s = 1
    for k in range(1, n + 1):
        s *= k
    return s


## 2. Exponential Function
def exponential(x, tol=1e-14, max_iterations=1000):
    """Compute e**x using a Taylor series.

    We write x = x0 + z where x0 = int(round(x)). Then
    e**x = e**x0 * sum_{n=0} (z**n / n!). Stop when |s| < tol
    or after `max_iterations` terms.
    """
    if not isinstance(x, (int, float)):
        raise TypeError("exponential() argument must be a number") # updated to check type

    x0 = int(round(x)) # updated to use int(round(x))
    z = x - x0

    s = 1.0
    series_sum = s
    for n in range(1, max_iterations + 1):
        s = (z ** n) / factorial(n)
        series_sum += s
        if abs(s) < tol:
            break

    # scale by e^x0 
    return (e ** x0) * series_sum


## 3. Natural Logarithm Function
def natural_log(x, tol=1e-12, max_iterations=100):
    """Compute ln(x)

    Solve f(s) = e^s - x = 0 with Newton updates.
    """
    if not isinstance(x, (int, float)):
        raise TypeError("argument must be a number") # updated to check type
    if x <= 0:
        raise ValueError("input must be positive") # updated to check positivity

    # quick check
    if x == 1.0:
        return 0.0

    # initial guess
    s = x

    # Newton iteration: s_next = s - (e^s - x)/e^s = s - 1 + x * e^{-s}
    for k in range(max_iterations):
        s_next = s - 1 + x * exponential(-s)
        if abs(s_next - s) < tol:
            return s_next
        s = s_next

    raise RuntimeError("natural_log failed to converge")


## Part 2: Solve the logistic growth problem
"""
The population P(t) satisfies the logistic growth equation:
    dP/dt = r * P * (1 - P/K)
    where P(t) is the population at time t, with initial condition P(0) = P0, (0< P0 < K),
    r is the intrinsic growth rate, and K > 0 is the carrying capacity.
    
    The logistic equation admits a closed-form solution:
    P(t) = K / (1 + ((K - P0)/P0) * exponential(-r * t))
    
    1. Using the exact formula above, evaluate P(t) for a list of time values
    2. Implement the forward Euler method to approximate the solution of the logistic equation
    3. Find the time when the population hits K/2
    
    For numerical experiments, use r = 0.5, K = 100, PO = 10, and t ranges from 0 to 20
"""
## Build an absolute value function
def abs_val(x):
    """Return the absolute value of x."""
    if x < 0:
        return -x
    else:
        return x
    
## Build a square root function using Newton's method
def sqrt(x, init_guess=1.0, tol=1e-14, max_iterations=100):
    """Compute the square root of x using Newton's method."""
    if x < 0:
        raise ValueError("Cannot compute square root of negative number") # check for negative input
    
    s = init_guess
    for _ in range(max_iterations):
        s_new = 0.5 * (s + x / s)
        if abs_val(s_new - s) < tol:
            return s_new
        s = s_new
    raise RuntimeError("sqrt failed to converge")

### Quick test of abs_val and sqrt functions
# if __name__ == "__main__":
#      print("abs_val(-5) =", abs_val(-5))
#      print("sqrt(25) =", sqrt(25))

## Let us proceed to the main logistic growth problem
def population_exact(t, P0=10, r=0.5, K=100):
    """Compute the exact population P(t) using the closed-form solution. I will attempt this by breaking down the formula
    into smaller parts for clarity."""
    exponent = -r * t
    exp_value = exponential(exponent)
    denominator = 1 + ((K - P0) / P0) * exp_value
    P_t = K / denominator
    return P_t

## Quick test of population_exact function
# if __name__ == "__main__":
#     t_test = 10.0
#     print(f"Exact population at time t={t_test} is P(t)={population_exact(t_test)}")

## Finding the time when population hits K/2 using natural_log
def time_to_half_population(P0=10, r=0.5, K=100):
    """Compute the time t when population P(t) hits K/2 using the closed-form solution."""
    half_K = K / 2
    ratio = (K - P0) / P0
    ln_argument = ratio
    ln_value = natural_log(ln_argument)
    t_half = (1 / r) * ln_value
    return t_half

## Quick test of time_to_half_population function
# expected answer for defaults: t=4.394449154672438


## Implementing the forward Euler method
"""
Plan:
1. Find h, which is 20/N. Ensure N is a range so that we can loop through it and it automatically adjusts h.
    So we do for N in N_values:
        h = 20/N ...
2. We have to check that the absolute value of numerical solution - exact solution is less than tol. 
    Tol is a constant multiplied by h. Such as abs_val(numerical - exact) < C1 * h ...
3. We then need to find the time when the absolute value of numerical solution - K/2 is less than tol.
    Such as abs_val(numerical - K/2) < C2 * h ...
4. Set t to go from 0 to 20
"""
def population_euler(t_start=0.0, t_end=20.0, N=1000, P0=10, r=0.5, K=100):
    """Simple forward Euler integrator returning time and population arrays.

    Returns (t_values, P_values).
    """
    if N <= 0:
        raise ValueError("N must be a positive integer")

    h = (t_end - t_start) / N
    t_vals = [t_start]
    P_vals = [P0]

    t = t_start
    P = P0
    for _ in range(N):
        P = P + h * r * P * (1 - P / K)
        t = t + h
        t_vals.append(t)
        P_vals.append(P)

    return t_vals, P_vals


def find_crossing_time(t_vals, P_vals, target):
    """Find first time when P crosses `target` (in this case, K/2) by linear interpolation between steps.

    Returns None if crossing not found.
    """
    for i in range(1, len(P_vals)):
        P0 = P_vals[i - 1]
        P1 = P_vals[i]
        if (P0 - target) == 0:
            return t_vals[i - 1]
        if (P0 - target) * (P1 - target) < 0: # crossing detected
            # linear interpolation
            t0 = t_vals[i - 1]
            t1 = t_vals[i]
            frac = (target - P0) / (P1 - P0)
            return t0 + frac * (t1 - t0)
    return None


def run_logistic_comparison(t_start=0.0, t_end=20.0, N=1000, P0=10, r=0.5, K=100, tol_factor=1.0):
    """Run a comparison between the exact solution and forward Euler.

    - Computes Euler solution with N steps
    - Evaluates exact solution at the Euler time points
    - Computes max absolute error and checks against `tol`
    - Finds crossing times for P=K/2 for both solutions
    Returns a dict of results.
    """
    t_euler, P_euler = population_euler(t_start, t_end, N, P0, r, K)

    # tolerance scaled with step size h
    h = (t_end - t_start) / N
    tol = tol_factor * h

    # evaluate exact solution at Euler time points for error computation
    P_exact_at_euler = []
    for i in range(len(t_euler)):
        ti = t_euler[i]
        val = population_exact(ti, P0, r, K)
        P_exact_at_euler.append(val)

    max_error = 0.0
    for i in range(len(P_euler)):
        pn = P_euler[i]
        pe = P_exact_at_euler[i]
        err = abs_val(pn - pe)
        if err > max_error:
            max_error = err

    within_tol = max_error <= tol

    t_half_exact = time_to_half_population(P0, r, K)
    t_half_euler = find_crossing_time(t_euler, P_euler, K / 2)

    results = {
        "t_euler": t_euler,
        "P_euler": P_euler,
        "P_exact_at_euler": P_exact_at_euler,
        "max_error": max_error,
        "within_tol": within_tol,
        "t_half_exact": t_half_exact,
        "t_half_euler": t_half_euler,
        "params": {"t_start": t_start, "t_end": t_end, "N": N, "P0": P0, "r": r, "K": K, "tol": tol},
    }

    return results

# Running the pipeline with a test case
if __name__ == "__main__":
    res = run_logistic_comparison(N=10000000, tol_factor=100) # Tolerance factor increased for demonstration
    print(f"Parameters: {res['params']}")
    print(f"Max absolute error between Euler and exact: {res['max_error']}")
    print(f"Within tolerance: {res['within_tol']}")
    print(f"Analytic time to K/2: {res['t_half_exact']}")
    print(f"Euler estimated time to K/2: {res['t_half_euler']}")
