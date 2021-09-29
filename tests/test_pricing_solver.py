from math import floor, log10, lcm
import decimal
from collections import defaultdict

from pygcgopt import GCGModel, quicksum, PricingSolver, GCG_PRICINGSTATUS

from ortools.algorithms import pywrapknapsack_solver

from test_cpmp import build_model

import pytest

import os


class PyKnapsackSolver(PricingSolver):
    def solve(self, pricingprob, probnr, dualsolconv):
        vars = {var.name: var for var in pricingprob.getVars()}
        if pricingprob.getNBinVars() + pricingprob.getNIntVars() < len(vars):
            return {"status": GCG_PRICINGSTATUS.NOTAPPLICABLE}
        for var_name in vars:
            if pricingprob.isNegative(vars[var_name].getLbLocal()):
                return {"status": GCG_PRICINGSTATUS.NOTAPPLICABLE}
        conss = pricingprob.getConss()
        if len(conss) != 1:
            return {"status": GCG_PRICINGSTATUS.NOTAPPLICABLE}
        
        cons = conss[0]

        conshdlr_name = cons.getConstraintHandlerName()
        if conshdlr_name == "linear":
            if not pricingprob.isIntegral(pricingprob.getRhs(cons)) or not pricingprob.isInfinity(-pricingprob.getLhs(cons)):
                return {"status": GCG_PRICINGSTATUS.NOTAPPLICABLE}

            consvals = pricingprob.getValsLinear(cons)
            if not all(pricingprob.isIntegral(v) for v in consvals.values()):
                return {"status": GCG_PRICINGSTATUS.NOTAPPLICABLE}
            consvals = {k: floor(v) for k, v in consvals.items()}
            
            capacity = floor(pricingprob.getRhs(cons))

            prelcapacity = capacity
            inferbounds = False
            for var_name, consval in consvals.items():
                if pricingprob.isInfinity(vars[var_name].getUbLocal()):
                    inferbounds = True
                
                if consval < 0:
                    if pricingprob.isInfinity(vars[var_name].getUbLocal()):
                        return {"status": GCG_PRICINGSTATUS.NOTAPPLICABLE}
                    prelcapacity -= floor(consval * vars[var_name].getUbLocal())
            
            ubs = {}
            for var_name, consval in consvals.items():
                if inferbounds and pricingprob.isInfinity(vars[var_name].getUbLocal()):
                    ubs[var_name] = floor(abs(prelcapacity / consval))
                else:
                    ubs[var_name] = vars[var_name].getUbLocal()
        elif conshdlr_name == "knapsack":
            raise NotImplementedError("knapsack constraints are not implemented")
        else:
            return {"status": GCG_PRICINGSTATUS.NOTAPPLICABLE}

        item_var_map = {}
        profits = []
        weights = []
        for var_name, consval in consvals.items():
            items_per_var = int(ubs[var_name] - vars[var_name].getLbLocal() + 0.5)
            for i in range(items_per_var):
                item_var_map[len(profits)] = var_name
                profits.append(-vars[var_name].getObj())
            if pricingprob.isGE(vars[var_name].getLbLocal(), 1) and not pricingprob.isEQ(vars[var_name].getUbLocal(), 0.0):
                capacity -= floor(vars[var_name].getLbLocal() * consval)
        # print(f"{profits = }")
        
        # Adjust profits to be integer, multiply by lcm
        # calculating true precision using log is too slow, just assume 10^-6
        # prec = int(log10(round(1/pricingprob.epsilon())))
        prec = 6
        profits_lcm = lcm(*[decimal.Context(prec=prec).create_decimal_from_float(p).as_integer_ratio()[1] for p in profits])
        profits_lcm = min(profits_lcm, 1000000)
        profits = [int(p*profits_lcm) for p in profits]

        # Ugly solution, might not be the best in all cases
        # profits = [round(p*1000) for p in profits]

        for idx in range(len(item_var_map)):
            var_name = item_var_map[idx]
            if consvals[var_name] > 0:
                weights.append(consvals[var_name])
            else:
                capacity -= consvals[var_name]
                weights.append(-consvals[var_name])
                profits[idx] *= -1

        solitems = []
        if capacity < 0:
            return {"status": GCG_PRICINGSTATUS.INFEASIBLE}
        elif capacity == 0:
            solitems = []
        else:
            knapsackSolver = pywrapknapsack_solver.KnapsackSolver(
                pywrapknapsack_solver.KnapsackSolver.KNAPSACK_DYNAMIC_PROGRAMMING_SOLVER, "KnapsackExample"
            )
            knapsackSolver.Init(profits, [weights], [int(capacity)])
            knapsackSolver.Solve()
            solitems = [idx for idx in range(len(item_var_map)) if knapsackSolver.BestSolutionContains(idx)]

        solvals = defaultdict(int)
        for idx in range(len(item_var_map)):
            var_name = item_var_map[idx]
            assert(consvals[var_name] >= 0 or not vars[var_name].isNegated())
            if idx in solitems:
                if consvals[var_name] >= 0 and not vars[var_name].isNegated():
                    solvals[var_name] += 1.0
            else:
                if consvals[var_name] < 0 or vars[var_name].isNegated():
                    solvals[var_name] += 1.0
        
        for var_name, var in vars.items():
            if pricingprob.isGE(var.getLbLocal(), 1):
                solvals[var_name] += floor(var.getLbLocal())

        col = pricingprob.createGcgCol(probnr, [vars[var_name] for var_name in solvals], [solvals[var_name] for var_name in solvals], False, pricingprob.infinity())
        self.model.getMasterProb().addCol(col)

        sol_val = sum([solval * vars[var_name].getObj() for var_name, solval in solvals.items()])

        return {"lowerbound": sol_val, "status": GCG_PRICINGSTATUS.OPTIMAL}

    def solveHeuristic(self, pricingprob, probnr, dualsolconv):
        return self.solve(pricingprob, probnr, dualsolconv)


@pytest.mark.parametrize("lp_file,dec_file", [
    ("instances_bpp/N1C1W4_M.BPP.lp", "instances_bpp/N1C1W4_M.BPP.dec"),
    ("instances_bpp/N1C2W2_O.BPP.lp", "instances_bpp/N1C2W2_O.BPP.dec")
])
def test_pypricer_fast(lp_file, dec_file):
    dirname = os.path.dirname(__file__)
    lp_file = os.path.join(dirname, lp_file)
    dec_file = os.path.join(dirname, dec_file)

    m = GCGModel()
    m.readProblem(lp_file)
    m.readProblem(dec_file)

    m.optimize()

    assert(m.getStatus() == "optimal")

    gcg_pricer_sol_obj_val = m.getSolObjVal(m.getBestSol())

    m = GCGModel()
    m.readProblem(lp_file)
    m.readProblem(dec_file)

    for p in m.listPricingSolvers():
        m.setPricingSolverEnabled(p, False)

    proxyKnapsackSolver = PyKnapsackSolver()
    m.includePricingSolver(proxyKnapsackSolver, "pyknapsack", "Python ortools knapsack pricing solver", 300, False, True)

    m.optimize()

    assert(m.getStatus() == "optimal")

    py_pricer_sol_obj_val = m.getSolObjVal(m.getBestSol())

    assert(gcg_pricer_sol_obj_val == py_pricer_sol_obj_val)
