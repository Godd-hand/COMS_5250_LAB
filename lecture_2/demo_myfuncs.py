# Python Lab Assignment

# constant e
e = 2.7182818284590451


# 1. Factorial
def factorial(n):
    """Compute n! for a non-negative integer n.

    Assumes `n` is an integer >= 0.
    """
    s = 1
    for k in range(1, n + 1):
        s *= k
    return s


# 2. Exponential Function
def exponential(x, tol=1e-12, max_iterations=1000):
    """Compute e**x using a Taylor series.

    We write x = x0 + z where x0 = int(round(x)). Then
    e**x = e**x0 * sum_{n=0} (z**n / n!). Stop when |s| < tol
    or after `max_iterations` terms.
    """
    if not isinstance(x, (int, float)):
        raise TypeError("exponential() argument must be a number")

    x0 = int(round(x))
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


# 3. Natural Logarithm Function
def natural_log(x, tol=1e-12, max_iterations=100):
    """Compute ln(x)

    Solve f(s) = e^s - x = 0 with Newton updates.
    """
    if not isinstance(x, (int, float)):
        raise TypeError("argument must be a number")
    if x <= 0:
        raise ValueError("input must be positive")

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


if __name__ == "__main__":
    # demonstration 

    print("factorial(5) =", factorial(5))
    print("expected factorial(5) = 120")

    for val in (1.0, 2.5, -1.5, 10.0):
        try:
            print(f"exponential({val}) =", exponential(val))
            print(f"expected e**{val} =", (e ** val))
        except Exception as err:
            print("exponential error for", val, "->", err)

    for val in (e, 10.0, 0.5):
        try:
            print(f"natural_log({val}) =", natural_log(val))
        except Exception as err:
            print("natural_log error for", val, "->", err)

