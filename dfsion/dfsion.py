import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import matplotlib.animation as animation
plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10)

with open("zadanie1.dat", 'r') as file:
    header = list(map(float, file.readline().strip().split()))
plansza_raw = np.loadtxt("zadanie1.dat", skiprows=1)

def make_circle(x_0, y_0, r):
    def circle(x, y):
        return (x - x_0)**2 + (y - y_0)**2 < r**2
    return circle


theta = np.linspace(0, 2*np.pi, 100)

obszar_x, obszar_y, obszar_r = header[4], header[5], header[6]
obszarx = obszar_x + obszar_r*np.cos(theta)
obszary = obszar_y + obszar_r*np.sin(theta)

znikanie_x, znikanie_y, znikanie_r = header[7], header[8], header[9]
znikaniex = znikanie_x + znikanie_r*np.cos(theta)
znikaniey = znikanie_y + znikanie_r*np.sin(theta)

fig, ax = plt.subplots(figsize = (4,4))
ax.plot(plansza_raw[:,0], plansza_raw[:,1], 'o', color='k', alpha=1, markersize=1.5)
ax.plot(obszarx, obszary, '--', color='blue', label='obszar dyfuzji')
ax.plot(znikaniex, znikaniey, '--', color='y', label='obszar znikania')
ax.plot(header[0], header[1], 'o', color='r', label='punkt kreacji', alpha = 0.7)
ax.set_aspect('equal')
#ax.legend(frameon = False)

plt.tight_layout()
plt.title(r"N$_{max}$ = " + f'{int(round(header[3]))}, N = {int(round(header[12]))}', fontsize = 10)
plt.savefig(r"omega" + f"{int(header[2]/header[11])}" + r"Ra" + f"{header[9]}" + f"N{int(round(header[12]))}.pdf", format='pdf')
plt.show()

with open("zadanie1_wyniki.dat", 'r') as file:
    header = list(map(float, file.readline().strip().split()))
plansza_raw = np.loadtxt("zadanie1_wyniki.dat", skiprows=1).transpose()

av_Dxx = np.mean(plansza_raw[1])
av_Dyy = np.mean(plansza_raw[2])
av_Dxy = np.mean(plansza_raw[3])

av_Dxx2 = np.mean(plansza_raw[1]**2)
av_Dyy2 = np.mean(plansza_raw[2]**2)
av_Dxy2 = np.mean(plansza_raw[3]**2)

sigma_Dxx = np.sqrt((av_Dxx2 - av_Dxx**2)/len(plansza_raw[1]))
sigma_Dyy = np.sqrt((av_Dyy2 - av_Dyy**2)/len(plansza_raw[2]))
sigma_Dxy = np.sqrt((av_Dxy2 - av_Dxy**2)/len(plansza_raw[3]))

linsp = np.linspace(0, header[12]*header[11], len(plansza_raw[0]))

tits = [
    r"D$_{xx}$",
    r"D$_{yy}$",
    r"D$_{xy}$"
]

fig1, ax1 = plt.subplots(figsize=(4, 4))
ax1.plot(linsp, header[10]*np.ones_like(linsp), '--', color='k')
ax1.plot(linsp, np.zeros_like(linsp), '--', color='k')
ax1.plot(linsp, plansza_raw[1], label=tits[0], linewidth=0.7)
ax1.plot(linsp, plansza_raw[2], label=tits[1], linewidth=0.7)
ax1.plot(linsp, plansza_raw[3], label=tits[2], linewidth=0.7)

ax1.set_title(r"N$_{max}$ = " + f'{int(round(header[3]))}')
ax1.set_xlabel('t')
ax1.set_ylabel('$D$')
ax1.set_xlim(0, max(linsp))
ax1.legend(frameon=True, framealpha=1, edgecolor='k', facecolor='white', fancybox=False)
ax1.grid(True)
plt.tight_layout()
#plt.savefig("zad2.pdf", format = 'pdf')
plt.show()

tit2 = r"$\omega$" + f" = {int(header[2]/header[11])}" + r", R$_{a}$ = " + f"{header[9]}"
tit2_short = r"omega" + f"{int(header[2]/header[11])}" + r"Ra" + f"{header[9]}"

fig2, ax2 = plt.subplots(figsize=(4, 4))
ax2.plot(linsp, plansza_raw[0], color='k', linewidth=0.7)
ax2.set_title('liczba aktywnych cząstek')
ax2.set_xlabel('t')
ax2.set_ylabel(r'n')
ax2.set_title(tit2)
ax2.grid(True)
plt.tight_layout()
plt.savefig(f"{tit2_short}.pdf", format = 'pdf')
plt.show()

aaa = [
    r"D$_{xx}$ = " + f"{av_Dxx:.6f}" + r", $\sigma_{D_{xx}}$ = " + f"{sigma_Dxx:.6f}",
    r"D$_{yy}$ = " + f"{av_Dyy:.6f}" + r", $\sigma_{D_{yy}}$ = " + f"{sigma_Dyy:.6f}",
    r"D$_{xy}$ = " + f"{av_Dxy:.6f}" + r", $\sigma_{D_{xy}}$ = " + f"{sigma_Dxy:.6f}"
]

for tit in aaa:
    print(tit)