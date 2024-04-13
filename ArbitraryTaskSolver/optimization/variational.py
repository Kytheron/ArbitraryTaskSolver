def VTS_max_final_super_extra(Integralus, Boundaries):
    def mydisplay(expr):
        display(Latex(latex(expr)))
    t = Integralus.args[1][0]
    horizon = list(Integralus.args[1][1:])
    L = Integralus.args[0]
    variables = list(L.atoms(AppliedUndef))
    mydisplay(variables)
    mydisplay(L)
    mydisplay(Boundaries)

    def diff_degree_order(L):
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

    def get_lagrange_equations(L, variables):
        d_order = diff_degree_order(L)
        Lagrange_equations = []
        for i in range(len(variables)):
            LagrangeEq = 0
            for k in range(  max(list(d_order[i].values())[0] + 1,2)):
                a = Add( (-1)**k  *  (L.diff(variables[i].diff(t,k))).diff(t,k), LagrangeEq)
                LagrangeEq = a
            Lagrange_equations.append(Eq(LagrangeEq,0))
        return Lagrange_equations

    def sanitized_diffeqs(Lagrange):
        sanitized_diffeqs  = []
        for d in Lagrange:
            for v in variables:
                subs_dict = {k: Dummy('y') for k in d.find(Derivative)}
                d_ = d.subs(subs_dict)
                d_ = factor(d_.lhs - d_.rhs, v)
                wild = Wild('wild')
                if (wild * v).matches(d_):
                    d = simplify(Eq(d.lhs / v, d.rhs / v))
            sanitized_diffeqs.append(d)
        return sanitized_diffeqs

    def boundary_for_numeric(Boundaries):
        boundary_for_numeric = []
        for b in Boundaries:
          for i in range(len(b)):
            boundary_for_numeric.append(b[i])
        return boundary_for_numeric

    def find_consts_with_derivatives(Solution, B, variables):    
        const_eq = []
        Diff_Solution = [Eq(Solution[v].lhs.diff(t), Solution[v].rhs.diff(t)) for v in range(len(variables))]
        for v in range(len(variables)):
          for i in range(len(B[v])):
              if list(diff_degree_order(B[v][i])[0].values())[0] == 0:
                const_eq.append( Eq(B[v][i].rhs, Solution[v].rhs.subs(t,B[v][i].lhs.args[0])) )
              elif list(diff_degree_order(B[0][i])[0].values())[0] == 1:
                const_eq.append( Eq(B[v][i].rhs, Diff_Solution[v].rhs.subs(t,B[v][i].lhs.args[2][0])) )
        return solve(const_eq)

    def FinalSolution(Solution, const_sol, variables):
        FinalSolution = []
        for i in range(len(variables)):
            FinalSolution.append(Solution[i].subs(const_sol))
        return FinalSolution
    
    def fixed_solution(FinalSol,variables,B):
        FixedSolution = []
        for i in range(len(FinalSol)):
          for j in range(len(variables)):
            if FinalSol[i].rhs.subs(t,B[j][0].args[0].args[0]) == B[j][0].rhs:
              FixedSolution.append( Eq(B[j][0].lhs.subs(B[j][0].args[0].args[0],t),  FinalSol[i].rhs  ) )
        mydisplay(FixedSolution)
        return FixedSolution

    Lagrange = get_lagrange_equations(L,variables)
    mydisplay(Lagrange)
    boundary_for_numeric = boundary_for_numeric(Boundaries)
    counter = np.sum([1 for eq in Lagrange if eq.find(Derivative)])
    if counter == len(Lagrange): #diff system only
        try:
            if len(Lagrange) == 1:
                Solution = [dsolve(Lagrange[0])]
            else:
                Solution = dsolve(Lagrange, variables)
            mydisplay(Solution)
            const_sol = find_consts_with_derivatives(Solution, Boundaries, variables)
            mydisplay(const_sol)
            FinalSolution = FinalSolution(Solution, const_sol, variables)
            mydisplay(FinalSolution)
            FixedSolution = fixed_solution(FinalSolution, variables, Boundaries)
            lam_f = lambdify(t, FinalSolution[0].rhs)
            t = np.linspace(float(horizon[0]), float(horizon[1]), 100)
            sol = [ lam_f(t[i]) for i in range(len(t))]
            plt.plot(t, sol, color = 'blue', lw = 2 )
        except NotImplementedError:
              T, ret, vars_names = main(Lagrange, boundary_for_numeric,horizon[0], horizon[1])
              for i in range(len(ret)):
                plt.scatter(T, ret[vars_names[i]][1::], label=vars_names[i])

    elif counter == 0: #algebraic system
        if len(Lagrange) == 1:
            print(Lagrange[0])
            Solution = [solve(Lagrange[0])]
        else:
            Solution = solve(Lagrange, variables)
        print(Solution)
        const_sol = find_consts_with_derivatives(Solution, Boundaries, variables)
        FinalSolution = [i.subs(const_sol) for i in Solution]
        print(FinalSolution)

    T, ret, vars_names = main(Lagrange, boundary_for_numeric,float(horizon[0]), float(horizon[1]))
    for i in range(len(ret)):
      plt.scatter(T, ret[vars_names[i]][1::], label=vars_names[i])
    return Solution
g = VTS_max_final_super_extra(I, B)