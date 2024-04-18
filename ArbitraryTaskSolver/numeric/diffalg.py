import numpy as np
from sympy import Derivative, Function, Symbol, Eq, Integer, nsolve, Add, symbols, Wild
from sympy.core.function import AppliedUndef

t = symbols('t')
x = symbols('x', cls=Function)
def main(system, boundaries, tmin, tmax):
    def add_leading_nulls(s, N):
      return '0' * (N - len(s)) + s

    def naive_extrapolation(l, N):
      return l + [l[-1] for _ in range(N - len(l))]

    def diff_degree_order(L,variables):
        pairs_var_order = [{i.args[0] : i.args[1][1]} for i in list(L.find(Derivative))]
        variables_and_order = []
        for k in variables:
            diff_order = []
            for j in range(len(pairs_var_order)):
                if list(pairs_var_order[j].keys())[0] == k:
                    diff_order.append([list(pairs_var_order[j].values() )[0]])
            try:
                variables_and_order.append( {k : max(diff_order)[0] })
            except ValueError:
                variables_and_order.append( {k : 0 })
        return variables_and_order

    N = 50
    h = (tmax - tmin) / N
    finite_diff_dict = {
        # i == current state
        4: lambda y, h, i: (y.func(i+4) - 4 * y.func(i+3) + 6*y.func(i+2) - 4*y.func(i+1) + y.func(i))/ (h ** 4),
        3: lambda y, h, i: (y.func(i+3) - 3 * y.func(i+2) + 3*y.func(i+1) - y.func(i))/(h ** 3),
        2: lambda y, h, i: (y.func(i+2) - 2 * y.func(i+1) + y.func(i))/ (h ** 2),
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
    tau = {i.args[0] for i in diff_variables}.pop()

    # diff order total aldready done for diffeqs, now add variables
    diff_order_total = diff_order_total | {v: 0 for v in functions}
    variables = list(functions)

    # new arch
    i = Symbol('i')

    def gryaz(expr):
        return expr if expr in functions else expr.args[0]
    linearized_subs = {k: finite_diff_dict[v](gryaz(k), h, i) for k,v in diff_order_total.items()}
    linearized_subs |= {vremya: tmin + i*h}
    linearized_system = [eq.subs(linearized_subs) for eq in system]
    full_system = []
    setka_variables = set()
    left = 1
    right = N + 1 - diff_order + 2 ## в питоне надо делать +1 - (eq_degree +1)
    for k in range(left,right):
        for eq in linearized_system:
            linearized_i = eq.subs({i: k})#.subs({b.lhs.func(Integer((b.lhs.args[0] - tmin)/h)+1): b.rhs for b in boundaries})
            full_system.append(linearized_i)
            setka_variables.update(linearized_i.find(AppliedUndef))

    for bnd in boundaries:
        c1,c2,c3,c4 = symbols('c1 c2 c3 c4', cls=Wild)
        a,b,c,d = symbols('a b c d', cls=Wild)
        wild_const = [c1,c2,c3,c4]
        wild_third_boundaries = c1 + c2 + c3 + c4
        output_equation = 0
        for c in wild_const:
          if str(Eq(bnd.lhs.match(wild_third_boundaries)[c], 0)) == 'True':
            pass
          else:
            if diff_degree_order(Eq(bnd.lhs.match(wild_third_boundaries)[c], 0), [variables[0]] )[0][variables[0]] == 0:
              k0, dif0 = symbols('k0 dif0', cls = Wild)
              wild_zero_order = dif0 * k0
              ax_eq = Eq(bnd.lhs.match(wild_third_boundaries)[c], 0)
              q = Eq(ax_eq.lhs.match(wild_zero_order)[k0], 0)
              eq_to_add_0 = ax_eq.lhs.match(wild_zero_order)[dif0]*(q.lhs.func(Integer((q.lhs.args[0] - tmin)/h)+1))
              output_equation = Add(output_equation, eq_to_add_0)
            if diff_degree_order(Eq(bnd.lhs.match(wild_third_boundaries)[c], 0), [variables[0]] )[0][variables[0]] == 1:
              k1, dif1 = symbols('k1 dif1', cls = Wild)
              wild_first_order = dif1 * k1
              ax_eq = Eq(bnd.lhs.match(wild_third_boundaries)[c], 0)
              q = Eq(ax_eq.lhs.match(wild_first_order)[k1], 0)
              x = q.lhs.args[0].args[0]
              i = Integer((q.lhs.args[2][0] - tmin)/h + 1)
              eq_to_add_1 = ax_eq.lhs.match(wild_first_order)[dif1]*(x.subs(t, i+1) - x.subs(t, i) )/h
              output_equation = Add(output_equation, eq_to_add_1)
            if diff_degree_order(Eq(bnd.lhs.match(wild_third_boundaries)[c], 0), [variables[0]] )[0][variables[0]] == 2:
              k2, dif2 = symbols('k2 dif2', cls = Wild)
              wild_second_order = dif2 * k2
              ax_eq = Eq(bnd.lhs.match(wild_third_boundaries)[c], 0)
              q = Eq(ax_eq.lhs.match(wild_second_order)[k2], 0)
              x = q.lhs.args[0].args[0]
              i = Integer((q.lhs.args[2][0] - tmin)/h + 1)
              eq_to_add_2 = ax_eq.lhs.match(wild_second_order)[dif2]*(x.subs(t, i+2) - 2*x.subs(t, i+1) + x.subs(t, i)  )/h**2
              output_equation = Add(output_equation, eq_to_add_2)
            if diff_degree_order(Eq(bnd.lhs.match(wild_third_boundaries)[c], 0), [variables[0]] )[0][variables[0]] == 3:
              k3, dif3 = symbols('k3 dif3', cls = Wild)
              wild_third_order = dif3 * k3
              ax_eq = Eq(bnd.lhs.match(wild_third_boundaries)[c], 0)
              q = Eq(ax_eq.lhs.match(wild_third_order)[k3], 0)
              x = q.lhs.args[0].args[0]
              i = Integer((q.lhs.args[2][0] - tmin)/h + 1)
              eq_to_add_3 = ax_eq.lhs.match(wild_third_order)[dif3]*(x.subs(t, i+3) - 3*x.subs(t, i+2) + 3*x.subs(t, i+1) - x.subs(t, i)  )/h**3
              output_equation = Add(output_equation, eq_to_add_3)
        full_system.append(Eq(output_equation, bnd.rhs))

    sorted_setka_variables = sorted(setka_variables,key=lambda i: i.name + '_' + add_leading_nulls(i.args[0].__str__(), len(str(N))), reverse=False)
    vars_names = sorted(set([i.name for i in sorted_setka_variables]), reverse=False)
    sols = nsolve(full_system, sorted_setka_variables, [0 for _ in setka_variables])
    mapping = {k: v for k,v in zip(sorted_setka_variables, sols)}
    ret = {vars_names[i]: naive_extrapolation([v for k,v in mapping.items() if vars_names[i] in k.name], N+1) for i in range(len(vars_names))}
    T = np.linspace(tmin, tmax, N+1)

    return T, ret, vars_names