import numpy as np
import matplotlib.pyplot as plt

# Time array
t = np.linspace(0, 10, 1000)

# Frequencies (in THz or arbitrary units)
f_ir = 3.0      # frequency of IR radiation
f_dipole = 3.0  # same frequency => resonance case
# Try f_dipole = 2.5 to see off-resonance

# Waves
ir_wave = np.sin(2 * np.pi * f_ir * t)
dipole_wave = 0.8 * np.sin(2 * np.pi * f_dipole * t)

# Set dark background and light text to match LaTeX theme
plt.style.use('dark_background')
plt.rcParams.update({
    'figure.facecolor': '#282A36',
    'axes.facecolor': '#282A36',
    'axes.edgecolor': '#D9D9D9',
    'axes.labelcolor': '#D9D9D9',
    'text.color': '#D9D9D9',
    'xtick.color': '#D9D9D9',
    'ytick.color': '#D9D9D9',
    'grid.color': '#44475A',  # slightly lighter than bg
    'legend.facecolor': '#282A36',
    'legend.edgecolor': '#D9D9D9',
})

# Plot the waves
plt.figure(figsize=(10, 6))
plt.plot(t, ir_wave, label='IR Radiation (Electric Field)', color='#00BFFF', linewidth=2)
plt.plot(t, dipole_wave, label='Molecular Dipole Moment', color='#FF6347', linestyle='--', linewidth=2)

# Highlight overlap region (energy absorption)
plt.fill_between(t, ir_wave, dipole_wave, color='#50FA7B', alpha=0.3, label='Energy Transfer (Absorption)')

# Labels and styling
plt.title("IR Absorption — When IR Frequency Matches the Changing Dipole Moment", fontsize=14)
plt.xlabel("Time (a.u.)", fontsize=12)
plt.ylabel("Amplitude", fontsize=12)
plt.legend(fontsize=11, loc='upper right')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
