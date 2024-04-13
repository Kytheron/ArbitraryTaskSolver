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