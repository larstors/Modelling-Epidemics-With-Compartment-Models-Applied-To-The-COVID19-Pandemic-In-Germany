
import numpy as np
import scipy.integrate as sci
import time 
import threading
import sys
import multiprocessing as multiproc
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool 

# the function library
from projectlib import *
import h5py


# for timing
t3 = time.time()

# array for time
t = np.zeros(100000)

# matrices for different solves (with different parameters)
sol1 = np.array([np.zeros(100000) for i in range(4*38)])
sol2 = np.array([np.zeros(100000) for i in range(4*38)])
sol3 = np.array([np.zeros(100000) for i in range(4*38)])
sol4 = np.array([np.zeros(100000) for i in range(4*38)])
sol5 = np.array([np.zeros(100000) for i in range(4*38)])

# file with the numerical data
file = np.loadtxt("./NumData/alpha.txt")


# filling the arrays
k = 0
for row in file:
    t[k] = row[0]
    for i in range(4*38):
        sol1[i, k] = row[i+1]
        sol2[i, k] = row[i+1+38*4]
        sol3[i, k] = row[i+1+8*38]
        sol4[i, k] = row[i+1+12*38]
        sol5[i, k] = row[i+1+16*38]

    k += 1


def plotting(index):
    """Function for plotting the effect that alpha has on the differential equations

    Args:
        index (int): index of cell
    """
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(9,7), sharex=True, sharey=True)
    fig.suptitle(r"Graphical presentation the differential equation in cell %g $\beta = 1/14$, $p = 0.0264$, $t_0 = 10/24$" % (index + 1))

    axs[0, 0].set_title("S")
    axs[0, 0].set_xlabel("Time in days")
    axs[0, 0].set_ylabel("Part of population")
    axs[0, 0].plot(t, sol1[index, :], "b", label=r"$alpha = 0.15$")
    axs[0, 0].plot(t, sol2[index, :], "r", label=r"$alpha = 0.2$")
    axs[0, 0].plot(t, sol3[index, :], "g", label=r"$alpha = 0.24$")
    axs[0, 0].plot(t, sol4[index, :], "k", label=r"$alpha = 0.28$")
    axs[0, 0].plot(t, sol5[index, :], "m", label=r"$alpha = 0.4$")
    axs[0, 0].legend(loc="upper right")
    axs[0, 0].grid()


    axs[0, 1].set_title("I")
    axs[0, 1].set_xlabel("Time in days")
    axs[0, 1].set_ylabel("Part of population")
    axs[0, 1].plot(t, sol1[index + 38, :], "b", label=r"$alpha = 0.15$")
    axs[0, 1].plot(t, sol2[index + 38, :], "r", label=r"$alpha = 0.2$")
    axs[0, 1].plot(t, sol3[index + 38, :], "g", label=r"$alpha = 0.24$")
    axs[0, 1].plot(t, sol4[index + 38, :], "k", label=r"$alpha = 0.28$")
    axs[0, 1].plot(t, sol5[index + 38, :], "m", label=r"$alpha = 0.4$")
    axs[0, 1].legend(loc="upper right")
    axs[0, 1].grid()


    axs[1, 0].set_title("R")
    axs[1, 0].set_xlabel("Time in days")
    axs[1, 0].set_ylabel("Part of population")
    axs[1, 0].plot(t, sol1[index + 2*38, :], "b", label=r"$alpha = 0.15$")
    axs[1, 0].plot(t, sol2[index + 2*38, :], "r", label=r"$alpha = 0.2$")
    axs[1, 0].plot(t, sol3[index + 2*38, :], "g", label=r"$alpha = 0.24$")
    axs[1, 0].plot(t, sol4[index + 2*38, :], "k", label=r"$alpha = 0.28$")
    axs[1, 0].plot(t, sol5[index + 2*38, :], "m", label=r"$alpha = 0.4$")
    axs[1, 0].legend(loc="upper left")
    axs[1, 0].grid()


    axs[1, 1].set_title("D")
    axs[1, 1].set_xlabel("Time in days")
    axs[1, 1].set_ylabel("Part of population")
    axs[1, 1].plot(t, sol1[index + 3*38, :], "b", label=r"$alpha = 0.15$")
    axs[1, 1].plot(t, sol2[index + 3*38, :], "r", label=r"$alpha = 0.2$")
    axs[1, 1].plot(t, sol3[index + 3*38, :], "g", label=r"$alpha = 0.24$")
    axs[1, 1].plot(t, sol4[index + 3*38, :], "k", label=r"$alpha = 0.28$")
    axs[1, 1].plot(t, sol5[index + 3*38, :], "m", label=r"$alpha = 0.4$")
    axs[1, 1].legend(loc="upper right")
    axs[1, 1].grid()

    plt.savefig("Media/Plots/alpha/alphaplot_%g.pdf" % index, dpi=200)
    plt.close(fig)




if __name__ == "__main__":
    """for i in range(38):
        plotting(i)
    """
    
    pool = Pool(4)
    pool.map(plotting, range(38))
    
    """

    arg = [(i,) for i in range(38)]

    print(tuple(arg))




    pool = ThreadPool(2)
    pool.starmap(plotting, arg) 
    pool.close() 
    pool.join()
    """
    
t4 = time.time()
print("The program took ", t4- t3, "s")


"""
hf = h5py.File("./NumData/alpha.h5", "r")

print(list(hf.keys()))

n1 = hf.get("data")

n1 = np.array(n1)

print(n1)"""