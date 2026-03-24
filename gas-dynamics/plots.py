import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

loc = ".\\zadanie_gaz1\\zad1\\10_5\\"
data0 = np.loadtxt(loc + "hist_0.dat").T
data1 = np.loadtxt(loc + "hist_1.dat").T
data2 = np.loadtxt(loc + "hist_2.dat").T
data3 = np.loadtxt(loc + "hist_3.dat").T
data10 = np.loadtxt(loc + "hist_10.dat").T
data100 = np.loadtxt(loc + "hist_100.dat").T
data1000 = np.loadtxt(loc + "hist_1000.dat").T
data_fin = np.loadtxt(loc + "hist_fin.dat").T

data_list = [data0, data1, data2, data3, data10, data100, data1000, data_fin]
iteration_labels = ['0', '1', '2', '3', '10', '100', '1000', 'ostatnia']

n = len(data_list)
ncols = 2
nrows = (n + 1)//ncols

fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3*nrows), sharex=True)
axes = axes.flatten()

plot_idx = 0
for idx, (d, label) in enumerate(zip(data_list, iteration_labels)):
    if idx == 4:
        continue
    if plot_idx == n - 1:
        break 
    ax = axes[plot_idx]
    observed = d[1]
    expected = d[2]
    mask = expected > 0
    chi2 = np.sum((observed[mask] - expected[mask])**2 / expected[mask])
    ax.bar(d[0], d[1], width=(d[0][1] - d[0][0]), label='histogram', alpha=1, color='g')
    ax.plot(d[0], d[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
    ax.set_title(f"iteracja: {label}\n" + r"$\chi^2$" + f"={chi2:.2e}")
    ax.set_xlabel('v (m/s)')
    ax.set_ylabel('p-stwo')
    ax.legend(frameon = False)
    plot_idx += 1

ax = axes[n - 1]
ax.bar(data0[0], data0[1], width=(data0[0][1] - data0[0][0]), label='na początku', alpha=0.5, color='g')
ax.bar(data_fin[0], data_fin[1], width=(data_fin[0][1] - data_fin[0][0]), label='ostatnia iteracja', alpha=0.5, color='b')
ax.set_title("początek i koniec")
ax.set_xlabel('v (m/s)')
ax.set_ylabel('p-stwo')
ax.legend(frameon = False)

for j in range(n, len(axes)):
    fig.delaxes(axes[j])

fig.tight_layout()
plt.show()

loc = ".\\zadanie_gaz1\\zad4\\10_6\\"
data0 = np.loadtxt(loc + "hist_0.dat").T
data1 = np.loadtxt(loc + "hist_1.dat").T
data2 = np.loadtxt(loc + "hist_2.dat").T
data3 = np.loadtxt(loc + "hist_3.dat").T
data10 = np.loadtxt(loc + "hist_10.dat").T
data100 = np.loadtxt(loc + "hist_100.dat").T
data1000 = np.loadtxt(loc + "hist_1000.dat").T
data_fin = np.loadtxt(loc + "hist_fin.dat").T

data_list = [data0, data1, data2, data3, data10, data100, data1000, data_fin]

fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

axs[0].bar(data0[0], data0[1], width=(data0[0][1] - data0[0][0]), label='histogram', alpha=0.5, color='gray')
axs[0].plot(data0[0], data0[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
axs[0].set_title("iteracja: 1")
axs[0].set_xlabel('v (m/s)')
axs[0].set_ylabel('p-stwo')
axs[0].legend(frameon=False)

axs[1].bar(data_fin[0], data_fin[1], width=(data_fin[0][1] - data_fin[0][0]), label='histogram', alpha=0.5, color='gray')
axs[1].plot(data_fin[0], data_fin[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
axs[1].set_title("iteracja: 20000")
axs[1].set_xlabel('v (m/s)')
axs[1].set_ylabel('p-stwo')
axs[1].legend(frameon=False)

axs[2].bar(data100[0], data100[1], width=(data100[0][1] - data100[0][0]), label='histogram', alpha=0.5, color='gray')
axs[2].plot(data100[0], data100[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
axs[2].set_title("iteracja: 100")
axs[2].set_xlabel('v (m/s)')
axs[2].set_ylabel('p-stwo')
axs[2].legend(frameon=False)

axs[3].bar(data1000[0], data1000[1], width=(data1000[0][1] - data1000[0][0]), label='histogram', alpha=0.5, color='gray')
axs[3].plot(data1000[0], data1000[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
axs[3].set_title("iteracja: 1000")
axs[3].set_xlabel('v (m/s)')
axs[3].set_ylabel('p-stwo')
axs[3].legend(frameon=False)

plt.tight_layout()
plt.show()

fig2, axs2 = plt.subplots(2, 2, figsize=(12, 10))
axs2 = axs2.flatten()

axs2[0].bar(data0[0], data0[1], width=(data0[0][1] - data0[0][0]), label='histogram', alpha=0.5, color='gray')
axs2[0].plot(data0[0], data0[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
axs2[0].set_title("iteracja: 1")
axs2[0].set_xlabel('v (m/s)')
axs2[0].set_ylabel('p-stwo')
axs2[0].legend(frameon=False)

axs2[1].bar(data1[0], data1[1], width=(data1[0][1] - data1[0][0]), label='histogram', alpha=0.5, color='gray')
axs2[1].plot(data1[0], data1[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
axs2[1].set_title("iteracja: 1")
axs2[1].set_xlabel('v (m/s)')
axs2[1].set_ylabel('p-stwo')
axs2[1].legend(frameon=False)

axs2[2].bar(data2[0], data2[1], width=(data2[0][1] - data2[0][0]), label='histogram', alpha=0.5, color='gray')
axs2[2].plot(data2[0], data2[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
axs2[2].set_title("iteracja: 2")
axs2[2].set_xlabel('v (m/s)')
axs2[2].set_ylabel('p-stwo')
axs2[2].legend(frameon=False)

axs2[3].bar(data_fin[0], data_fin[1], width=(data_fin[0][1] - data_fin[0][0]), label='histogram', alpha=0.5, color='gray')
axs2[3].plot(data_fin[0], data_fin[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
axs2[3].set_title("iteracja: 20000")
axs2[3].set_xlabel('v (m/s)')
axs2[3].set_ylabel('p-stwo')
axs2[3].legend(frameon=False)

plt.tight_layout()
plt.show()

fig3, ax3 = plt.subplots(figsize=(8, 5))
ax3.bar(data0[0], data0[1], width=(data0[0][1] - data0[0][0]), label='pierwsza iteracja', alpha=0.5, color='k')
ax3.bar(data_fin[0], data_fin[1], width=(data_fin[0][1] - data_fin[0][0]), label='ostatnia iteracja', alpha=0.5, color='r')
ax3.plot(data0[0], data0[2], label='r. Maxwella-Boltzmanna', linestyle='-', color='k')
ax3.set_title("iteracja: początek i koniec")
ax3.set_xlabel('v (m/s)')
ax3.set_ylabel('p-stwo')
ax3.legend(frameon=False)
plt.show()

loc = ".\\zadanie_gaz1\\zad3\\10_5\\"
data0 = np.loadtxt(loc + "\\rv_0.dat").T
data1 = np.loadtxt(loc + "\\rv_1.dat").T
data2 = np.loadtxt(loc + "\\rv_2.dat").T
data3 = np.loadtxt(loc + "\\rv_3.dat").T
data10 = np.loadtxt(loc + "\\rv_10.dat").T
data100 = np.loadtxt(loc + "\\rv_100.dat").T
data1000 = np.loadtxt(loc + "\\rv_1000.dat").T
data_fin = np.loadtxt(loc + "\\rv_fin.dat").T

data_sets = [data0, data1, data10, data1000]
titles = ["iteracja: 0", "iteracja: 2", "iteracja: 10", "iteracja: ostatnia"]

fig4, ax4 = plt.subplots(2, 2, figsize=(10, 10))
ax4 = ax4.flatten()

for i, (data, title) in enumerate(zip(data_sets, titles)):
    ax_left = ax4[i]
    norm = (data[2] - data[2].min()) / (data[2].max() - data[2].min())
    blackness = 0.5 + 0.5 * norm  # 0.5 (gray) to 1 (black)
    colors = [(1 - b, 1 - b, 1 - b) for b in blackness]  # RGB for grayscale

    ax_left.set_xticks([])
    ax_left.set_yticks([])
    ax_left.scatter(data[0], data[1], color=colors, label='data[4]')
    ax_left.set_title(title)
    ax_left.tick_params(axis='y', labelcolor='tab:blue')
    ax_left.set_xlim(np.min(data1000[0]), np.max(data1000[0]))
    ax_left.set_ylim(np.min(data1000[1]), np.max(data1000[1]))

plt.tight_layout()
plt.show()

loc = ".\\zadanie_gaz1\\zad4\\10_6"
datazero = np.loadtxt(loc + "\\wyniki\\nptv_0.dat").T
data100 = np.loadtxt(loc + "\\wyniki\\nptv_100.dat").T
data1000 = np.loadtxt(loc + "\\wyniki\\nptv_10000.dat").T
data_fin = np.loadtxt(loc + "\\wyniki\\nptv_19900.dat").T
datas = [datazero, data100, data1000, data_fin]
tits = ['0', '100', '10000', '19900']

fig, axs = plt.subplots(1, 2, figsize=(8, 4), sharex=True)
for idx, data0 in enumerate(datas):
    
    axs[0].plot(data0[0], data0[3], marker='o', label =f"it.: {tits[idx]}")
    axs[1].plot(data0[0], data0[2], marker='o', label = f"it.: {tits[idx]}")

axs[0].legend(frameon = False)
axs[1].legend(frameon = False)
axs[0].set_xlabel('x')
axs[0].set_ylabel('T')
axs[0].set_xlim(0, 1)
axs[0].set_ylim(bottom=0)
axs[1].set_xlabel('x')
axs[1].set_ylabel('p')
axs[1].set_xlim(0, 1)
axs[1].set_ylim(bottom=0)
plt.tight_layout()
plt.show()