# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:21:22 2025

@author: aayan
"""

import math as m
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.special as sp
from scipy.optimize import curve_fit

# Define the HLN function
def HLN(B, alpha, l_phi):
    e = 1.602e-19
    pi = np.pi
    h = 6.62607015e-34
    hbar = h / (2 * pi)

    B = np.abs(B) + 1e-12  # avoid division by zero
    B_phi = hbar / (4 * e * l_phi**2)
    
    psi_term = sp.digamma(0.5 + B_phi / B)
    log_term = np.log(B_phi / B)
    
    return alpha * e**2 / (2 * pi**2 * hbar) * (psi_term - log_term) * h / e**2


# List of temperatures and filenames
temperatures = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])  # Numeric temperature values
temp_labels = ['2K', '3K', '4K', '5K', '6K', '7K', '8K', '9K', '10K']
exp_colors = ["red", "green", "blue", "magenta", "cyan", "purple", "brown", "pink", "gray"]
fit_colors = ["black", "darkred", "darkgreen", "darkblue", "navy", "indigo", "maroon", "darkorange", "gold"]

alpha_values = []  # Store fitted alpha values
l_phi_values = []  # Store fitted l_phi values

plt.figure(figsize=(12, 8))

# Loop through each temperature
for temp, label, exp_color, fit_color in zip(temperatures, temp_labels, exp_colors, fit_colors):
    filename = f"{label}.CSV"
    df = pd.read_csv(filename)

    # For 2K, the MR is in the 4th column (index 3), else 3rd column (index 2)
    col_idx = 3 if label == '2K' else 2
    B = df.iloc[:, 1].values * 1e-4  # convert from Gauss to Tesla
    R = df.iloc[:, col_idx].values


    delta_MC =  (1 / R - np.min(1/R))
    safe_B = np.abs(B)

    # Curve fitting with only l_phi bounded
    popt, _ = curve_fit(HLN, safe_B, delta_MC, p0=[10, 10e-9], maxfev=5000000)
    
    alpha_fit, l_phi_fit = popt
    alpha_values.append(alpha_fit)
    l_phi_values.append(l_phi_fit * 1e9)  # Convert L_phi from meters to nm
    
    delta_MC_fit = HLN(safe_B, *popt)

    # Compute R² value (correctly squared differences)
    ss_res = np.sum((delta_MC - delta_MC_fit) ** 2)
    ss_tot = np.sum((delta_MC - np.mean(delta_MC)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    # Output fit parameters clearly
    print(f"{label}: α = {alpha_fit:.6f}, L_φ = {l_phi_fit:.2e} m, R² = {r_squared:.6f}")

    # Plot data and fits correctly
    plt.scatter(B, delta_MC, label=f"{label} data", color=exp_color, s=10)
    plt.plot(B, delta_MC_fit, linewidth=1.8, color=fit_color, linestyle="-", label=f"{label} HLN Fit")

plt.xlabel(r"$H$ (Tesla)", fontsize=14)
plt.ylabel(r"$\Delta\sigma$ ($e^2/h$)", fontsize=14)
plt.legend(fontsize=10, ncol=2)
plt.tick_params(direction='in', length=5, width=1) 
plt.legend(fontsize='small', markerscale=0.6, bbox_to_anchor=(1.2, 1), loc='upper left')
plt.title("HLN Fit at Different Temperatures", fontsize=16)
plt.tight_layout()
plt.show()

fig, ax1 = plt.subplots(figsize=(8, 6))


ax1.plot(temperatures, alpha_values, marker='o', linestyle='-', color='b', label=r'$\alpha$')
ax1.set_xlabel('Temperature (K)', fontsize=14)
ax1.set_ylabel(r'$\alpha$', color='b', fontsize=14)
ax1.tick_params(axis='y', labelcolor='b')

ax2 = ax1.twinx()
ax2.plot(temperatures, l_phi_values, marker='s', linestyle='-', color='r', label=r'$L_{\phi}$ (nm)')
ax2.set_ylabel(r'$L_{\phi}$ (nm)', color='r', fontsize=14)
ax2.tick_params(axis='y', labelcolor='r')


ax1.legend(loc='upper left', fontsize=12)
ax2.legend(loc='upper right', fontsize=12)

plt.title(r'$\alpha$ and $L_{\phi}$ vs Temperature', fontsize=16)
plt.show()
