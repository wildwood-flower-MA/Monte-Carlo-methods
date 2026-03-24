import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import matplotlib.animation as animation
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)

name = str(sys.argv[1])
plansza_raw = np.loadtxt("lab3.dat")
size = np.shape(plansza_raw)

const = (np.max(plansza_raw[:,0]) - np.min(plansza_raw[:,0]))/(np.max(plansza_raw[:,1]) - np.min(plansza_raw[:,1]))


fig, ax = plt.subplots(figsize = (10,10))
ax.plot(plansza_raw[:,0], plansza_raw[:,1], 'o', color = 'k', alpha = 0.3)
ax.set_aspect('equal')
#ax.set_xlim((-1,1)); ax.set_ylim((-1,1))

plt.title(name, fontsize = 20)
plt.savefig(f"{name}.png")
plt.show()


