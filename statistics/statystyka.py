import numpy as np
import os
import matplotlib.pyplot as plt
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)

def f(x: float) -> float:
    return 0.8*(1.0 + x - x**3)

def F(x: float) -> float:
    return 0.8*x + 0.2*(2*(x**2) - x**4)

def F_1(x: float) -> float:
    return 0.8*x

def F_2(x: float) -> float:
    return 0.2*(2*(x**2) - x**4)

fig, axes = plt.subplots(1, 2, figsize=(15, 7))
axes[0].plot(np.linspace(0, 1, 100), f(np.linspace(0, 1, 100)), label=r'$f(x) = \frac{4}{5}(1+x-x^3)$', color='k')
axes[0].set_xlabel('x', fontsize=20)
axes[0].set_ylabel(r'$f(x)$', fontsize=20)
axes[0].grid(True)
axes[0].legend(fontsize=20)
axes[1].plot(np.linspace(0, 1, 100), F(np.linspace(0, 1, 100)), label=r'$F(x) = F_1(x)+F_2(x)$', color='k', linestyle='--')
axes[1].plot(np.linspace(0, 1, 100), F_1(np.linspace(0, 1, 100)), label=r'$F_1(x) = \frac{4}{5}x$', color='r', linestyle='--')
axes[1].plot(np.linspace(0, 1, 100), F_2(np.linspace(0, 1, 100)), label=r'$F_2(x) = \frac{1}{5}(2x^2-x^4)$', color='b', linestyle='--')
axes[1].set_xlabel('x', fontsize=20)
axes[1].grid(True)
axes[1].set_ylabel('$F(x)$', fontsize=20)
axes[1].legend(fontsize=20)
plt.savefig('f_F.png')
plt.tight_layout()
plt.show()

tabele = {}
dane = {}
linspaces = {}
etykiety = ["Rozkład złożony",
            r"Łańcuch Markowa, $\Delta = 0.5$",
            r"Łańcuch Markowa, $\Delta = 0.05$",
            "Metoda eliminacji"]
dane_raw = np.loadtxt("lab2.dat")

N = 10

for idx, etykieta in enumerate(etykiety):
    dane[etykieta] = dane_raw[idx,:-2]
    linspaces[etykieta] = np.linspace(dane_raw[idx,-2], dane_raw[idx,-1], 11)

for idx, etykieta in enumerate(etykiety):
    fig = plt.figure(figsize= (10,10))
    stairs_patch = plt.stairs(np.size(dane[etykieta])*dane[etykieta]/np.sum(dane[etykieta]), linspaces[etykieta],
            color = 'Gray', fill=True, label = "zliczenia", alpha = 0.5)
    plt.plot(linspaces[etykieta], f(linspaces[etykieta]),
        linestyle='-', color = 'k', label = "f-cja gęstości p-stwa")
    stairs_patch.set_edgecolor('darkgray')
    stairs_patch.set_linewidth(2)
    #plt.xlabel('X', fontsize = 20)
    plt.legend(fontsize = 20)
    #plt.ylabel("N", fontsize = 20)
    plt.ylim((0.8,1.2))
    plt.title(etykieta, fontsize = 20)
    plt.savefig(f"lab2_{idx}.png")
    plt.show()

chain = np.loadtxt("lab2_chain.dat")
fig, axes = plt.subplots(1, 2, figsize=(15, 7))


axes[0].plot(chain[:, 0], color='k', label=r'$\Delta = 0.5$')
axes[0].legend(fontsize=20)
axes[0].set_xlabel('numer losowania', fontsize=20, color='k')
axes[0].set_title(r'Łańcuch Markowa ($\Delta = 0.5$)', fontsize=20, color='k')
axes[0].tick_params(axis='both', labelsize=20)

axes[1].plot(chain[:, 1], color='r', label=r'$\Delta = 0.05$')
axes[1].legend(fontsize=20)
axes[1].set_xlabel('numer losowania', fontsize=20, color='k')
axes[1].set_title(r'Łańcuch Markowa ($\Delta = 0.05$)', fontsize=20, color='k')
axes[1].tick_params(axis='both', labelsize=20)

plt.tight_layout()
plt.show()