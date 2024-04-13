def main(system, boundaries, tmin, tmax):
    def add_leading_nulls(s, N):
      return '0' * (N - len(s)) + s

    def naive_extrapolation(l, N):
      return l + [l[-1] for _ in range(N - len(l))]

    N = 50
    h = (tmax - tmin) / N
    finite_diff_dict = {
        # i == current state
        4: lambda y, h, i: (y.func(i+2) - 4 * y.func(i+1) + 6*y.func(i) - 4*y.func(i-1) + y.func(i-2))/ (h ** 4),
        3: lambda y, h, i: (y.func(i+2) - 2 * y.func(i+1) + 2*y.func(i-1) - y.func(i-2))/ (2*h ** 3),
        2: lambda y, h, i: (y.func(i+1) - 2 * y.func(i) + y.func(i-1))/ (h ** 2),
        1: lambda y, h, i: (y.func(i+1) - y.func(i)) / h,
        0: lambda y, h, i: y.func(i)
    }
    # kraevaya zadacha?
    diff_variables = set()
    functions = set()
    diff_order_total = {}
    for eq in system:
        deriv_degree_dict = {i: i.args[1][1]  for i in [*eq.find(Derivative)]}
        diff_order = max([*deriv_degree_dict.values()]) if len(deriv_degree_dict) > 0 else 0
        diff_variable = [i.args[0] for i in deriv_degree_dict.keys()]
        diff_variables.update(diff_variable)
        diff_order_total = diff_order_total | deriv_degree_dict
        functions.update(eq.find(Function))

    vremya = {i.args[0] for i in diff_variables}.pop()

    variables = functions - diff_variables
    # diff order total aldready done for diffeqs, now add variables
    diff_order_total = diff_order_total | {v: 0 for v in functions}

    # new arch
    i = Symbol('i')

    def gryaz(expr):
        return expr if expr in functions else expr.args[0]
    linearized_subs = {k: finite_diff_dict[v](gryaz(k), h, i) for k,v in diff_order_total.items()}
    linearized_subs |= {vremya: tmin + i*h}
    linearized_system = [eq.subs(linearized_subs) for eq in system]
    #print(linearized_system)
    full_system = []
    setka_variables = set()
    for k in range(1,N+1):
        for eq in linearized_system:
            linearized_i = eq.subs({i: k})#.subs({b.lhs.func(Integer((b.lhs.args[0] - tmin)/h)+1): b.rhs for b in boundaries})
            full_system.append(linearized_i)
            setka_variables.update(linearized_i.find(AppliedUndef))

    full_system.extend([Eq(b.lhs.func(Integer((b.lhs.args[0] - tmin)/h)+1), b.rhs) for b in boundaries])
    sorted_setka_variables = sorted(setka_variables,key=lambda i: i.name + '_' + add_leading_nulls(i.args[0].__str__(), len(str(N))), reverse=False)
    vars_names = sorted(set([i.name for i in sorted_setka_variables]), reverse=False)
    sols = nsolve(full_system, sorted_setka_variables, [0 for _ in setka_variables])
    mapping = {k: v for k,v in zip(sorted_setka_variables, sols)}
    ret = {vars_names[i]: naive_extrapolation([v for k,v in mapping.items() if vars_names[i] in k.name], N+1) for i in range(len(vars_names))}
    T = np.linspace(tmin, tmax, N+1)
    return T, ret, vars_names