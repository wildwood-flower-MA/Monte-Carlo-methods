import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import matplotlib.animation as animation
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)

with open("temp.dat", "r") as temp:
    raw_data = np.loadtxt("temp.dat")
    labels = ["całka po A","całka po A","całka po B","całka po B"]
    fig, axes = plt.subplots(1, 2, figsize=(10,5), squeeze=False)

    ax1 = axes[0, 0]
    for i in [0, 2]:
        x = raw_data[i*3]
        y = raw_data[i*3 + 1]
        yerr = raw_data[i*3 + 2]
        ax1.errorbar(x, y, yerr=yerr, fmt='o', label=labels[i], capsize=5)
    ax1.set_xscale('log')
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'$10^{{{int(np.log10(val))}}}$' for val in x], fontsize=15)
    ax1.set_xlabel('n', fontsize=20)
    ax1.set_ylabel('wsp. pow.', fontsize=20)
    ax1.legend(fontsize=20)
    ax1.set_title(r'$x_a = R_b+0.5R_a$', fontsize=20)
    ax1.grid(True, which="both", linestyle='--', linewidth=0.5)

    ax2 = axes[0, 1]
    for i in [1, 3]:
        x = raw_data[i*3]
        y = raw_data[i*3 + 1]
        yerr = raw_data[i*3 + 2]
        ax2.errorbar(x, y, yerr=yerr, fmt='o', label=labels[i], capsize=5)
    ax2.set_xscale('log')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'$10^{{{int(np.log10(val))}}}$' for val in x], fontsize=15)
    ax2.set_xlabel('n', fontsize=20)
    ax2.set_ylabel('wsp. pow.', fontsize=20)
    ax2.legend(fontsize=20)
    ax2.set_title(r'$x_a = 0$', fontsize=20)
    ax2.grid(True, which="both", linestyle='--', linewidth=0.5)

    plt.tight_layout()

    plt.show()
