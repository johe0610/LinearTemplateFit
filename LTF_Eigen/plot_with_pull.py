# plot_with_pull.py
import numpy as np
import matplotlib.pyplot as plt

# Example data (replace with your actual arrays)

y = np.array([
    172.892, 172.865, 172.870, 172.834, 172.735,
    173.317, 172.689, 172.549, 173.111, 172.684,
    175.701, 173.137, 171.618, 172.961, 171.725,
    171.108, 157.861
])

yerr = np.array([2.12677,
                1.63042,
                1.73072,
                2.11719,
                3.3496,
                11.4183,
                1.72564,
                4.76383,
                22.9485,
                62.3693,
                36.1656,
                15.0441,
                17.398,
                18.1851,
                17.4037,
                12.4368,
                374.681
                ])
# Bin indices
x = np.arange(len(y))

# Compute pull: (BinContent - 172.75) / BinError
pull = (y - 172.75) / yerr

# Create figure with two subplots (main + pull panel)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
ax1.set_ylim(152.75, 192.75)
ax2.set_ylim(-2, 2)
# --- Top panel: Data with error bars and horizontal line ---
ax1.errorbar(x, y, yerr=yerr, fmt='o', color='blue', label='Data', markersize=4, capsize=3)
ax1.axhline(y=172.75, color='red', linestyle='--', linewidth=2, label='172.75')
ax1.set_ylabel('Value')
ax1.set_title('Data with Pull Panel')
ax1.legend()
ax1.grid(True, alpha=0.3)

# --- Bottom panel: Pull values ---
ax2.axhline(y=0, color='black', linestyle='-', linewidth=1)  # Zero line
ax2.errorbar(x, pull, fmt='o', color='purple', markersize=4, capsize=3)
ax2.set_ylabel('Pull')
ax2.set_xlabel('Bin Index')
ax2.grid(True, alpha=0.3)

# Adjust layout to prevent overlap
plt.tight_layout()

# Save as high-quality PDF (vector format)
plt.savefig("plot_with_pull.pdf", bbox_inches='tight', dpi=300)

# Show the plot (optional)
#plt.show()
