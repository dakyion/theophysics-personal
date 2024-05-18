## Written by Abdennour Harche
## May 2024

# Run this script to graph the variation of the Ground Energy with respect to the pairing Constant G.
# A comparison graph between classical variational method and VQE


import numpy as np
import matplotlib.pyplot as plt

nstates = input('Number of states: ')
nparticles= input('Number of particles: ')
g1 = np.linspace(0.1,1,10)
g2 = np.linspace(0.1,1,32)
y1 = np.fromfile(f'./BCSEnergiesS/Sbcs{nstates}s{nparticles}p{g1[0]}to{g1[-1]}g{len(g1)}s.dat', sep='\n')
#y2 = np.fromfile(f'./fermionthesis/ExactEnergies/exact{int(nstates)//2}s{nparticles}p{gg[0]}to{gg[-1]}g{len(gg)}s.dat', sep='\n')
y3 = np.fromfile(f'./EE/ee{int(nstates)-2}s{nparticles}p{g2[0]}to{g2[-1]}g{len(g2)}s.dat', sep='\n')

plt.plot(g1,y1,'+r', g2,y3,'-g')
plt.legend(['VQE-BCS','Classic-BCS'])
plt.xlabel(r'$G/\hbar\omega$')
plt.ylabel(r'$E_0/\hbar\omega$')
plt.show()