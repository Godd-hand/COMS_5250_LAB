import numpy as np
import matplotlib.pyplot as plt
import os

def PlotPoly():
    filename = 'poly.data'
    
    # check if file exists
    if not os.path.exists(filename):
        print(f"Error: {filename} not found")
        return

    # Read data
    try:
        data = np.loadtxt(filename)
        x = data[:, 0]
        y = data[:, 1]
    except Exception as e:
        print(f"Error reading data: {e}")
        return

    # Plotting configuration
    plt.rc("font", size=14)
    plt.figure(1, figsize=(8, 6))
    
    plt.plot(x, y, linestyle="dashed", linewidth=2, marker="o", color="blue", markersize=8)
    
    plt.xlim(-1.0, 1.0)
    plt.grid(True)
    
    plt.xlabel("x-axis", size=16)
    plt.ylabel("y-axis", size=16)
    plt.title("Chebyshev Polynomial Plot", size=18)
    
    output_file = 'chebyshev_plot.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.show()

if __name__ == "__main__":
    PlotPoly()