import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import matplotlib.animation as animation
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)

with open("temp.dat", "r") as temp:
    raw_data = np.loadtxt("temp.dat")

kolory = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']
fig, ax = plt.subplots(figsize = (10,10))
ax.set_aspect('equal')
for i in range(len(raw_data)//2):
    ax.plot(raw_data[2*i], raw_data[2*i+1], 'o', color = kolory[i], alpha = 0.3)
plt.title(r"Koła", fontsize = 20)
plt.show()