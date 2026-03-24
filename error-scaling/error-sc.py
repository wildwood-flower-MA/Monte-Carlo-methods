import numpy as np
import os
import matplotlib.pyplot as plt
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)

tabele = {}
dane = {}
etykiety = ['X', 'X2', 'err_X', 'var_num', 'var_teo', 'err_var', 'n']
etykiety_plus = ['01', '05', '09']
dane_raw = np.loadtxt("lab1.dat")

for idx, etykieta_plus in enumerate(etykiety_plus):
    dane[etykieta_plus] = dane_raw[idx*6:(idx+1)*6,:]

    for idxx, etykieta in enumerate(etykiety):
        tabele[etykieta+etykieta_plus] = dane[etykieta_plus][:,idxx]

    for idx, etykieta in enumerate(['err_X', 'err_var']):
        fig = plt.figure(figsize= (10,10))
        for idx, etykieta_plus in enumerate(etykiety_plus):
            plt.plot(tabele['n'+etykieta_plus],tabele[etykieta+etykieta_plus], linestyle='-', marker='o',
                    label = r"p = "+etykieta_plus[0]+r"."+etykieta_plus[1])
            plt.yscale('log'); plt.xscale('log')
            plt.xlabel('N', fontsize = 20)
        plt.legend(fontsize = 20)
        plt.ylabel(etykieta, fontsize = 20)
        plt.show()
        fig.savefig(etykieta+etykieta_plus)

fig = plt.figure(figsize= (10,5))
arr = np.linspace(0.005,0.995,100)
def err_var(beta: float ,arr: np.ndarray) -> np.ndarray:
    return np.abs(1-(beta/(arr-arr*arr)))

for b in [0.07,.1,.25,.49,0.81]:
    plt.plot(arr, err_var(beta = b, arr = arr), linestyle='-',
            label = r"$\beta$ = "+f'{b}')
    plt.ylim((0,8));
plt.xlabel('p', fontsize = 20)
plt.legend(fontsize = 20)
plt.ylabel('err_var', fontsize = 20)
plt.show()
fig.savefig('a')