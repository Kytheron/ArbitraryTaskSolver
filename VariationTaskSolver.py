class Updateable(object):
    def update(self, new):
        for key, value in new.__dict__.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
@dataclass
class OutputTemplate(Updateable):
    formula: list() #list of Eq, add typing here
    values: list()
    t_list: list()
    range: list()

def list_reinit(obj):
    if isinstance(obj, list):
        return obj
    else:
        return [obj]

def get_euler_lagrange_equations(L,variables):
    EulerLagrangeEqs = []
    for var in variables:
        EulerLagrangeEqs.append( Eq(L.diff(var), (L.diff(var.diff(t))).diff(t))  )
    return EulerLagrangeEqs

def constants_solution(Solution, boundary_conditions):
    ret = []
    solution_transform = {s.lhs.func: s.rhs for s in Solution}
    for b in boundary_conditions:
        subs_dict = {fnc: solution_transform[fnc.func].subs(t, fnc.args[0]) for fnc in b.find(AppliedUndef)}
        ret.append(b.subs(subs_dict))
    const_sol = solve(ret, dict = True)
    return const_sol[0]

def VTS_max_final(Integralus, Boundaries, forced_numeric=True, check_sync=True):
    t = Integralus.args[1][0]
    horizon = list(Integralus.args[1][1:])
    L = Integralus.args[0]
    variables = list(L.atoms(AppliedUndef))
    free_symbols = L.free_symbols - {t}
    Lagrange = get_euler_lagrange_equations(L,variables)
    counter = np.sum([1 for eq in Lagrange if eq.find(Derivative)])

    config_dict = {}
    if free_symbols:
        for v in free_symbols:
            config_dict[v] = input(f'Введите числовое значение для {v}:')

    final_numeric_ret = OutputTemplate(formula=None, values=None, t_list=None, range=horizon)
    if forced_numeric:
        T, ret, vars_names = main([i.subs(config_dict) for i in Lagrange], Boundaries, float(horizon[0]), float(horizon[1]))
        final_numeric_ret = OutputTemplate(formula=None, values=ret, t_list=T, range=horizon)


    if counter == len(Lagrange): #diff system only
        Solution = dsolve(Lagrange, variables)

        const_sol = constants_solution(Solution, Boundaries)
        FinalSolution = [i.subs(const_sol) for i in Solution]
        final_symbolic_ret = OutputTemplate(formula=FinalSolution, values=None, range=horizon, t_list=None)

    elif counter == 0: #algebraic system
        Solution = solve(Lagrange, variables)
        const_sol = constants_solution(Solution, Boundaries)
        FinalSolution = [i.subs(const_sol) for i in Solution]
        final_symbolic_ret = OutputTemplate(formula=FinalSolution, values=None, range=horizon, t_list=None)

    ret = OutputTemplate(formula=None, values=None, t_list=None, range=horizon)
    ret.update(final_symbolic_ret)


    if final_numeric_ret.values:
        ret.update(final_numeric_ret)
    else:
        applied_func_vals = {}
        for symb_ret in ret.formula:
            symb_ret_func = lambdify(t, symb_ret.rhs.subs(config_dict))
            symb_ret_values = symb_ret_func(T)
            applied_func_vals[symb_ret.lhs.func.name] = symb_ret_values
        final_apply_symbolic_ret = OutputTemplate(formula=None, values=applied_func_vals, t_list=None, range=horizon)
        ret.update(final_apply_symbolic_ret)


    if check_sync and forced_numeric:
        cum_mse_error = 0
        T = ret.t_list
        for symb_ret in ret.formula:
            symb_ret_func = lambdify(t, symb_ret.rhs.subs(config_dict))
            symb_ret_values = symb_ret_func(T)
            cum_mse_error += np.square(np.subtract(symb_ret_values, ret.values[symb_ret.lhs.func.name][:len(T)])).mean()

        print('Cummulative MSE error of symbolic and numeric outputs:', cum_mse_error)

    return ret

def visualize(obj: OutputTemplate):
    ngraphs = len(obj.values)
    t = obj.t_list
    if ngraphs > 1:
        fig, axs = plt.subplots(ngraphs, 1)
        for i, its in enumerate(obj.values.items()):
            axs[i].plot(t, its[1][:len(t)])
            axs[i].set_title(its[0])
        fig.show()
    else:
        its = [*obj.values.items()][0]
        fig = plt.plot(t, its[1][:len(t)])
        plt.title(its[0])
        plt.show()

ret = VTS_max_final(I, [*chain(*B)], forced_numeric=True)
visualize(ret)