
from matplotlib import pyplot as plt

from ArbitraryTaskSolver.numeric.diff import shooting_method
from ArbitraryTaskSolver.numeric.diff import numeric_cauchy_dsolve
from sympy import *
import numpy as np

from ArbitraryTaskSolver.numeric.diffalg import general_numerical_solver
from ArbitraryTaskSolver.optimization.variational import VTS


def shooting_method_test():
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

#if __name__ == '__main__':
#    shooting_method_test()

def numeric_dsolve_test():
    t = symbols('t')
    x = symbols('x', cls=Function)
    t_sol, x_sol = numeric_dsolve(0, 50, 500, Eq(x(t).diff(t, 2)  + sin(x(t)) + 75*x(t)**3, 0), [0, 1])
    plt.plot(t_sol, x_sol)
    plt.show()

#if __name__ == '__main__':
#    numeric_dsolve_test()

def main_test():
    t = symbols('t')
    x = symbols('x', cls=Function)
    # T, ret, vars_names = general_numerical_solver( [Eq( x(t).diff(t,2), 0  )], [ Eq( 2*(x(t).diff(t).subs(t,0)+ x(0)) , 1), Eq(x(0), 0) ] , 0, 1)
    # T, ret, vars_names = general_numerical_solver( [Eq( x(t).diff(t,3), 0  )], [ Eq( x(0) + 2*x(t).diff(t).subs(t,1), 4), Eq(x(1), 1), Eq(x(t).diff(t,2).subs(t,1), 2) ] , 0, 1)
    # T, ret, vars_names = general_numerical_solver( [Eq( x(t).diff(t,4), 0  )], [ Eq(x(0), 0), Eq(x(t).diff(t).subs(t,0), 0), Eq(x(t).diff(t,2).subs(t,1) + x(1), 8/3), Eq(x(t).diff(t,3).subs(t,0), 1) ] , 0, 1)
    # T, ret, vars_names = general_numerical_solver( [Eq( x(t).diff(t) - t*sin(t), 0  )], [ Eq(x(0), 0)] , 0, 1)
    T, ret, vars_names = general_numerical_solver([Eq(x(t).diff(t) - sin(t) ** 10 * cos(t), 0)], [Eq(x(0), 0)], 0, 5)
    for i in range(len(ret)):
        plt.scatter(T, ret[vars_names[i]][:len(T)], label=vars_names[i])
    time = np.linspace(0, 5, 500)
    plt.plot(time, np.sin(time) ** 11 / 11, color='blue')
    # plt.plot(time, np.sin(time) - time * np.cos(time) , color = 'blue')
    # plt.plot(time, time**3 / 6 + time**2 / 2, color = 'blue')
    # plt.plot(time, time**2, color = 'blue')
    plt.show()

#if __name__ == '__main__':
#    main_test()


def VTS_test():
    t = symbols('t')
    x = symbols('x', cls=Function)
    L = -2 * x(t) * t + (x(t).diff(t, 2)) ** 2
    I = Integral(L, (t, 0, 1))
    B = [[Eq(x(list(I.args[1])[1]), 0), Eq(x(list(I.args[1])[2]), 1 / 120), Eq((x(t).diff(t)).subs(t, list(I.args[1])[1]), 0), Eq((x(t).diff(t)).subs(t, list(I.args[1])[2]), 1 / 12)]]
    VTS(I, B)
    plt.show()

if __name__ == '__main__':
    VTS_test()