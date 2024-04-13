
from matplotlib import pyplot as plt

from ArbitraryTaskSolver.numeric.diff import shooting_method
from sympy import *
import numpy as np

def sosiska():
    tmin = 0
    tmax = 1
    t = symbols('t')
    x = symbols('x', cls=Function)
    N = 500
    alpha0 = 1.1899
    steps = 2
    def exact_sol():
        T = np.linspace(tmin,tmax,N)
        return T, np.pi * np.sin(T)/ np.sin(1)

    diffeq = Eq(x(t).diff(t, 2)  + x(t), 0)
    boundaries = [0, np.pi]

    t_sol, x_sol = shooting_method(tmin, tmax, N, steps, alpha0, diffeq, boundaries)

    plt.plot(t_sol, x_sol, color = 'green')
    plt.plot(*exact_sol(), color = 'red', lw = 2, ls = '--')
    plt.show()

if __name__ == '__main__':
    sosiska()