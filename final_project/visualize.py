import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def render_animation(global_nx=1152, global_ny=1152, size=4, num_steps=200000, interval=2000):
    local_nx = global_nx // size
    
    # Pre-calculate the meshgrid for the 3D surface plot
    x = np.arange(0, global_ny)
    y = np.arange(0, global_nx)
    X, Y = np.meshgrid(x, y)
    
    for step in range(0, num_steps + 1, interval):
        global_grid = np.zeros((global_nx, global_ny))
        
        # Stitch the distributed domains back together
        for rank in range(size):
            filename = f"output_rank{rank}_step{step:06d}.dat"
            if os.path.exists(filename):
                data = np.fromfile(filename, dtype=np.float64)
                data = data.reshape((local_nx, global_ny))
                global_grid[rank*local_nx:(rank+1)*local_nx, :] = data
                
        # Set up the 3D figure
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Render the surface. 
        # rstride and cstride downsample the polygons to prevent memory crashes on large grids.
        surf = ax.plot_surface(X, Y, global_grid, cmap='inferno', 
                               linewidth=0, antialiased=False, rstride=10, cstride=10)
        
        # Lock the Z-axis to the absolute min and max temperatures to keep the animation scale static
        ax.set_zlim(20, 5000)
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Temperature (°C)')
        plt.title(f'3D Heat Diffusion - Timestep {step:06d}')
        
        # Add a color bar
        fig.colorbar(surf, shrink=0.5, aspect=5, label='Temperature (°C)')
        
        # Save the frame with the 6-digit padding
        out_filename = f'frame_{step:06d}.png'
        plt.savefig(out_filename, dpi=150)
        plt.close()
        print(f"Rendered {out_filename}")

if __name__ == "__main__":
    # Ensure dimensions match the new divisible grid
    global_nx = 1152
    global_ny = 1152
    num_steps = 200000
    interval = 500
    
    # Auto-detect the number of ranks used in the last simulation run using 6-digit padding
    import glob
    rank_files = glob.glob("output_rank*_step000000.dat")
    size = len(rank_files)
    
    if size == 0:
        print("No output files found. Run the simulation first.")
    else:
        print(f"Detected {size} MPI ranks from output files. Starting 3D render...")
        render_animation(global_nx=global_nx, global_ny=global_ny, size=size, num_steps=num_steps, interval=interval)