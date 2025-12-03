cdef class PricingSolver:
    """Base class of the Pricing Solver Plugin"""

    cdef public Model model
    cdef public str solvername

    def freeSolver(self):
        pass

    def initSolver(self):
        pass

    def exitSolver(self):
        pass

    def initSolution(self):
        pass

    def exitSolution(self):
        pass

    def updateSolver(self, pricingprob, probnr, varobjschanged, varbndschanged, consschanged):
        pass

    def solve(self, pricingprob, probnr, dualsolconv):
        return {}

    def solveHeuristic(self, pricingprob, probnr, dualsolconv):
        return {}


cdef PricingSolver get_py_pricing_solver(GCG_SOLVER* pricingSolver) noexcept with gil:
    cdef GCG_SOLVERDATA* solverdata
    solverdata = GCGsolverGetData(pricingSolver)
    py_pricing_solver = <PricingSolver>solverdata
    return py_pricing_solver

cdef SCIP_RETCODE PyPricingSolverFree (GCG* gcg, GCG_SOLVER* solver) noexcept with gil:
    py_pricing_solver = get_py_pricing_solver(solver)
    py_pricing_solver.freeSolver()
    Py_DECREF(py_pricing_solver)
    return SCIP_OKAY

cdef SCIP_RETCODE PyPricingSolverInit (GCG* gcg, GCG_SOLVER* solver) noexcept with gil:
    py_pricing_solver = get_py_pricing_solver(solver)
    py_pricing_solver.initSolver()
    return SCIP_OKAY

cdef SCIP_RETCODE PyPricingSolverExit (GCG* gcg, GCG_SOLVER* solver) noexcept with gil:
    py_pricing_solver = get_py_pricing_solver(solver)
    py_pricing_solver.exitSolver()
    return SCIP_OKAY

cdef SCIP_RETCODE PyPricingSolverInitSol (GCG* gcg, GCG_SOLVER* solver) noexcept with gil:
    py_pricing_solver = get_py_pricing_solver(solver)
    py_pricing_solver.initSolution()
    return SCIP_OKAY

cdef SCIP_RETCODE PyPricingSolverExitSol (GCG* gcg, GCG_SOLVER* solver) noexcept with gil:
    py_pricing_solver = get_py_pricing_solver(solver)
    py_pricing_solver.exitSolution()
    return SCIP_OKAY

cdef SCIP_RETCODE PyPricingSolverUpdate (GCG* gcg, SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Bool varobjschanged, SCIP_Bool varbndschanged, SCIP_Bool consschanged) noexcept with gil:
    py_pricing_solver = get_py_pricing_solver(solver)
    py_pricingprob = GCGPricingModel.create(pricingprob)
    py_pricing_solver.updateSolver(py_pricingprob, probnr, varobjschanged, varbndschanged, consschanged)
    return SCIP_OKAY

cdef SCIP_RETCODE PyPricingSolverSolve (GCG* gcg, SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Real dualsolconv, SCIP_Real* lowerbound, GCG_PRICINGSTATUS* status) noexcept with gil:
    py_pricing_solver = get_py_pricing_solver(solver)
    py_pricingprob = GCGPricingModel.create(pricingprob)
    result_dict = py_pricing_solver.solve(py_pricingprob, probnr, dualsolconv)
    lowerbound[0] = result_dict.get("lowerbound", 0)
    status[0] = result_dict.get("status", <GCG_PRICINGSTATUS>status[0])
    return SCIP_OKAY

cdef SCIP_RETCODE PyPricingSolverSolveHeur (GCG* gcg, SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Real dualsolconv, SCIP_Real* lowerbound, GCG_PRICINGSTATUS* status) noexcept with gil:
    py_pricing_solver = get_py_pricing_solver(solver)
    py_pricingprob = GCGPricingModel.create(pricingprob)
    result_dict = py_pricing_solver.solveHeuristic(py_pricingprob, probnr, dualsolconv)
    lowerbound[0] = result_dict.get("lowerbound", 0)
    status[0] = result_dict.get("status", <GCG_PRICINGSTATUS>status[0])
    return SCIP_OKAY
