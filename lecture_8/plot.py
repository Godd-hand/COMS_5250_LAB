"""
Python code for plotting the results of lab assignment 8 (lab8.c).
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_results(filename):
    """
    Plot results for x and exp(x) as created in lab8.c.
    """
    # Load data from file
    data = np.loadtxt(filename)
    x_values = data[:, 0]
    exp_values = data[:, 1]

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, exp_values, label='exp(x)', color='blue', marker='o')
    plt.title('Plot of exp(x) vs x')
    plt.xlabel('x')
    plt.ylabel('exp(x)')
    plt.grid()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    plot_results('exp_data.txt')
    if len(sys.argv) > 1:
        plot_results(sys.argv[1])
    else:
        plot_results('exp_data.txt')