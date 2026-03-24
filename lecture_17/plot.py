import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    if not os.path.exists("original.dat") or not os.path.exists("blurred.dat"):
        print("Error: Image data files absent")
        return

    print("Loading image data in Python...")
    orig = np.loadtxt("original.dat")
    blur = np.loadtxt("blurred.dat")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # display original image
    ax1.imshow(orig, cmap='gray', vmin=0, vmax=255)
    ax1.set_title(f"Original Random Image\n({orig.shape[0]}x{orig.shape[1]})")
    ax1.axis('off')

    # display blurred image
    ax2.imshow(blur, cmap='gray', vmin=0, vmax=255)
    ax2.set_title(f"Blurred Image (3x3 Average)\n({blur.shape[0]}x{blur.shape[1]})")
    ax2.axis('off')

    plt.tight_layout()
    
    output_filename = "blur_comparison.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Visualization saved to {output_filename}")
    
    plt.show()

if __name__ == "__main__":
    main()