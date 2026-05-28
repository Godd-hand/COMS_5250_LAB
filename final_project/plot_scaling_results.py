import matplotlib.pyplot as plt
import numpy as np

# Nova timing data
cores = np.array([1, 2, 4, 8, 16, 32, 64])
times = np.array([439.15, 195.09, 103.81, 58.14, 30.42, 22.80, 27.34])

# Calculate speedup (T1 / Tp)
speedup = times[0] / times

# Apply a dark theme to match the 3D thermal diffusion animation
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(10, 7))

# Plot the ideal speedup reference line
ax.plot(cores, cores, linestyle='--', color='gray', label='Ideal Speedup', alpha=0.7)

# Plot the actual speedup curve
ax.plot(cores, speedup, marker='o', markersize=8, color='#d9f99d', linewidth=2.5, label='Actual Speedup')

# Highlight the peak efficiency at 32 cores
ax.plot(32, speedup[5], marker='o', markersize=14, markeredgecolor='#3b82f6', markerfacecolor='none', markeredgewidth=2)
ax.annotate(f'Peak Performance\n({speedup[5]:.1f}x at 32 Cores)', 
            xy=(32, speedup[5]), 
            xytext=(-45, 25),
            textcoords='offset points', 
            color='#3b82f6', 
            fontsize=11, 
            fontweight='bold',
            arrowprops=dict(arrowstyle='-', color='#3b82f6', alpha=0.5))

# Highlight the performance drop at 64 cores
ax.plot(64, speedup[6], marker='X', markersize=10, color='#ef4444')
ax.annotate('Communication\nOverhead', 
            xy=(64, speedup[6]), 
            xytext=(-45, -45),
            textcoords='offset points', 
            color='#ef4444', 
            fontsize=11, 
            fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#ef4444', alpha=0.8))

# Formatting the axes and title
ax.set_title('MPI Parallel Speedup on Nova Cluster', fontsize=16, fontweight='bold', pad=20)
ax.set_xlabel('Number of CPU Cores', fontsize=14, labelpad=10)
ax.set_ylabel('Speedup Factor', fontsize=14, labelpad=10)

# Set tick marks to match the tested core counts
ax.set_xticks(cores)
ax.set_ylim(0, 35)
ax.set_xlim(0, 68)

# Add a grid for readability
ax.grid(True, linestyle=':', alpha=0.3, color='white')
ax.legend(fontsize=12, loc='upper left', framealpha=0.2)

# Save the plot as a hi-res image
plt.tight_layout()
plt.savefig('scaling_plot_high_res.png', dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
print("Saved high-resolution plot to 'scaling_plot_high_res.png'")