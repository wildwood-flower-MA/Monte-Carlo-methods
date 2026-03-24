import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

potential = np.loadtxt("V_dokladnie.dat").T

plt.imshow(potential, cmap='Grays')
plt.colorbar(label='V')
#plt.title('nadrelaksacją')
plt.savefig("dokladnie.pdf")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.arange(potential.shape[1])
y = np.arange(potential.shape[0])
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, potential, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Potential')
plt.show()

from matplotlib import colors

potential1 = np.loadtxt("V_1.dat").T
potential2 = np.loadtxt("V_2.dat").T
potential3 = np.loadtxt("V_3.dat").T
potential1_diff = np.abs(potential1 - potential)
potential2_diff = np.abs(potential2 - potential)
potential3_diff = np.abs(potential3 - potential)
spotential1 = np.loadtxt("sigma_V_1.dat").T
spotential2 = np.loadtxt("sigma_V_2.dat").T
spotential3 = np.loadtxt("sigma_V_3.dat").T
s1 = np.loadtxt("S_1.dat").T
s2 = np.loadtxt("S_2.dat").T
s3 = np.loadtxt("S_3.dat").T

fig, ax = plt.subplots(4, 3, figsize = (12, 16))

column_labels = [r'N$_{ch.} = 100$, n$_{l} = 100$, B$ = 0$',
                r'N$_{ch.} = 100$, n$_{l} = 100$, B$ = 1$',
                r'N$_{ch.} = 300$, n$_{l} = 300$, B$ = 1$']
for i, label in enumerate(column_labels):
    ax[0, i].set_title(label)

row_labels = [r'V$_{MC}$', r'|V$_{MC}$ - V$_{NR}$|', r'$\sigma_V$', 'S']
for i, label in enumerate(row_labels):
    ax[i, 0].set_ylabel(label, rotation=90, size='large')

im00 = ax[0,0].imshow(potential1, cmap='Grays')
fig.colorbar(im00, ax=ax[0,0])
im01 = ax[0,1].imshow(potential2, cmap='Grays')
fig.colorbar(im01, ax=ax[0,1])
im02 = ax[0,2].imshow(potential3, cmap='Grays')
fig.colorbar(im02, ax=ax[0,2])

norm = colors.Normalize(vmin=0)
im10 = ax[1,0].imshow(potential1_diff, cmap='Grays', norm=norm)
fig.colorbar(im10, ax=ax[1,0])
im11 = ax[1,1].imshow(potential2_diff, cmap='Grays', norm=norm)
fig.colorbar(im11, ax=ax[1,1])
im12 = ax[1,2].imshow(potential3_diff, cmap='Grays', norm=norm)
fig.colorbar(im12, ax=ax[1,2])

im20 = ax[2,0].imshow(spotential1, cmap='Grays')
fig.colorbar(im20, ax=ax[2,0])
im21 = ax[2,1].imshow(spotential2, cmap='Grays')
fig.colorbar(im21, ax=ax[2,1])
im22 = ax[2,2].imshow(spotential3, cmap='Grays')
fig.colorbar(im22, ax=ax[2,2])

im30 = ax[3,0].imshow(s1, cmap='winter')
fig.colorbar(im30, ax=ax[3,0])
im31 = ax[3,1].imshow(s2, cmap='winter')
fig.colorbar(im31, ax=ax[3,1])
im32 = ax[3,2].imshow(s3, cmap='winter')
fig.colorbar(im32, ax=ax[3,2])

fig.tight_layout()
plt.savefig("wszystko.pdf")
plt.show()