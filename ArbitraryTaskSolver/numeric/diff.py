import numpy as np
import sympy
from sympy import Symbol, solve, Derivative


def numeric_dsolve(
        tmin, tmax, N, #сетка -- скаляры
        sympy_diffeq, #уравнение -- одно уравнение
        initial_boundaries # набор ограничение на задачу коши [в точке, первая произв, и тд] -- массив
    ):
    # шаг 1: валидация
    # уравнение: дифференц порядок уравнения >= 1
    deriv_degree_dict = {i: i.args[1][1]  for i in [*sympy_diffeq.find(Derivative)]}
    diff_order = max([*deriv_degree_dict.values()])
    main_variable = [*deriv_degree_dict.keys()][0].args[0] # костыль
    # print(main_variable)
    # cетка: tmax > tmin, N > порядок уравнения + 1
    assert tmax > tmin
    assert N > diff_order + 1
    # ограничения: длина массива огр == порядок уравнения
    assert len(initial_boundaries) == diff_order
    # шаг 2: решение
    h = (tmax - tmin)/N
    # преобразование диффура в разностную схему (от абстратных переменных)
    y = [Symbol(f'y_{i}') for i in range(diff_order + 1)]
    finite_diff_dict = {
        2: lambda y, h: (y[-1] - 2 * y[-2] + y[-3])/ (h ** 2),
        1: lambda y, h: (y[-2] - y[-3]) / h
    }
    linearized_subs = {k: finite_diff_dict[v](y, h) for k,v in deriv_degree_dict.items()} | {main_variable: y[-2]}
    ret = sympy_diffeq.subs(linearized_subs)
    y_next = solve(ret, y[-1])[0] # потому что одно решение для линейной сетки

    x_sol = [initial_boundaries[0], initial_boundaries[1] * h + initial_boundaries[0]] # реализовать для уравнений высших порядков
    for i in range(diff_order,N):
        x_sol.append(y_next.subs({y[0]: x_sol[i-2], y[1]: x_sol[i-1]}).subs({main_variable.args[0]: tmin + i*h}))

    t_sol = np.linspace(tmin, tmax, N)
    return t_sol, x_sol


def shooting_method(tmin, tmax, N, steps, alpha0, diffeq, boundaries):
    h = (tmax - tmin) / N
    t_sol, x_sol = numeric_dsolve(tmin, tmax, N, diffeq , [boundaries[0], alpha0])
    t_sol_h, x_sol_h = numeric_dsolve(tmin, tmax, N, diffeq, [boundaries[0], alpha0+h])
    alpha = [alpha0]
    for n in range(steps):
        alpha_n = alpha[-1] - ((x_sol[-1] - boundaries[1])*h / (x_sol_h[-1] - x_sol[-1]))
        alpha.append(alpha_n)
        t_sol, x_sol = numeric_dsolve(tmin, tmax, N, diffeq, [boundaries[0], alpha[-1]])
        t_sol_h, x_sol_h = numeric_dsolve(tmin, tmax, N, diffeq, [boundaries[0], alpha[-1]+h])

    return t_sol, x_sol

