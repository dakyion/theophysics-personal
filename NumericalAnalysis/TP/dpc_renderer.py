import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
file = np.loadtxt('./sp2000t10d8s9/step0.dat')
n = int(len(file)/2)
fig,ax = plt.subplots()
scat1 = ax.scatter(file[:n,0], file[:n,1], c='b',s=0.5)
scat2 = ax.scatter(file[n:,0], file[n:,1], c='r',s=0.5)
ax.set_xlim([-8e-7,2e-6])
ax.set_ylim([-8e-7,2e-6])
#plt.show()
# pv = file[:n]
# print(pv)
def update(frame):
    filename = "./sp2000t10d8s9/step{0}.dat".format(frame)
    file = np.loadtxt(filename)
    n = int(len(file)/2)
    ax.clear()
    ax.scatter(file[:n,0], file[:n,1], c='b',s=0.5)
    ax.scatter(file[n:,0], file[n:,1], c='r',s=0.5)
    ax.set_xlim([-8e-7,2e-6])
    ax.set_ylim([-8e-7,2e-6])
    del file
    return ax

ani = FuncAnimation(fig, update, frames=1000)
writer = animation.FFMpegWriter(fps=30,
                                metadata=dict(artist='Harche Abdennour'),
                                bitrate=1800)
ani.save('plasmasp2000t10d8s8.mp4', writer=writer)
# plt.show()
