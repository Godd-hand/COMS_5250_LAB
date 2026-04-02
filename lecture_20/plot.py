import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    if not os.path.exists("V_out.dat") or not os.path.exists("J_out.dat"):
        print("Error: Output data files not found") # some robustness in case files are missing
        return

    # 1. Plot V(x0)
    data_v = np.loadtxt("V_out.dat")
    x0 = data_v[:, 0]
    V_x0 = data_v[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(x0, V_x0, 'b-', linewidth=2)
    plt.title(r'1D Electrostatic Potential $V(x_0)$')
    plt.xlabel(r'$x_0$')
    plt.ylabel(r'$V(x_0)$')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig('V_potential.png', dpi=300, bbox_inches='tight')
    print("Saved Electrostatic Potential plot to V_potential.png")

    # 2. Plot J0(x)
    data_j = np.loadtxt("J_out.dat")
    x = data_j[:, 0]
    J_x = data_j[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(x, J_x, 'r-', linewidth=1.5)
    plt.axhline(0, color='black', linewidth=1)
    plt.title(r'Bessel Function $J_0(x)$')
    plt.xlabel('x')
    plt.ylabel(r'$J_0(x)$')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig('J_bessel.png', dpi=300, bbox_inches='tight')
    print("Saved Bessel Function plot to J_bessel.png")

    # Display plots
    plt.show()

if __name__ == "__main__":
    main()