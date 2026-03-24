import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import matplotlib.animation as animation
plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

data = np.loadtxt("histogram_wyniki_1_x1x2x3.txt").T
labels = [r"x$_1$", r"x$_2$", r"x$_3$"]
plt.figure(figsize = (6,5))
for i, x in enumerate(data[1:]):
    plt.plot(data[0], x, 'o', markersize=0.5, label=labels[i])
plt.xlabel("t")
plt.ylabel("x")
plt.legend(fancybox = False)
plt.title(r"P$_{max}$ = 1")
plt.grid()
plt.savefig("1.pdf")
plt.show()

data = np.loadtxt("histogram_wyniki_x1x2x3.txt").T
labels = [r"x$_1$", r"x$_2$", r"x$_3$"]
plt.figure(figsize = (6,5))
for i, x in enumerate(data[1:]):
    plt.plot(data[0], x, 'o', markersize=0.5, label=labels[i])
plt.xlabel("t")
plt.ylabel("x")
plt.legend(fancybox = False)
plt.title(r"P$_{max}$ = 5")
plt.grid()
plt.savefig("2.pdf")
plt.show()

data = np.loadtxt("histogram_wyniki_100.txt").T
t = data[0]
x3 = data[1]
sigma = data[2]

plt.errorbar(t, x3, yerr=sigma, fmt='o', markersize=1,
            ecolor='k', capsize=3, linestyle='-', color='blue')
plt.xlabel("t")
plt.ylabel(r"$\overline{x}_3$")
plt.title(r"P$_{max}$ = 100")
plt.grid()
plt.savefig("3.pdf")
plt.show()

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(t, data[2], color='k', label = r'$\sigma_{z}$')
ax.set_xlabel('t')
ax.set_ylabel(r'$\sigma_{3}$')
ax2 = ax.twinx()
ax2.plot(t, data[1], color='r', label = r'$\overline{x}_{3}$')
ax2.set_ylabel(r'$\overline{x}_{3}$')
fig.legend(loc='lower right', bbox_to_anchor=(0.8, 0.2))
fig.tight_layout()
plt.savefig("4.pdf")
plt.show()