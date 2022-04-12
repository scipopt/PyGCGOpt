from pyscipopt import Model, quicksum
import json

def build_model(n_locations, n_clusters, distances, demands, capacities):
    m = Model()

    m.printVersion()
    m.redirectOutput()

    m.setMinimize()

    x = {}
    y = {}

    for j in range(n_locations):
        y[j] = m.addVar(f"y_{j}", vtype="B")
        for i in range(n_locations):
            x[i, j] = m.addVar(f"x_{i}_{j}", vtype="B", obj=distances[min(i,j)][max(i,j)])

    # Hold different constraint types
    conss_assignment = []
    conss_capacity = []
    cons_pmedian = None

    # Create the assignment constraints
    for i in range(n_locations):
        conss_assignment.append(
            m.addCons(quicksum(x[i, j] for j in range(n_locations)) == 1)
        )

    # Create the capacity constraints
    for j in range(n_locations):
        conss_capacity.append(
            m.addCons(quicksum(demands[i] * x[i, j] for i in range(n_locations)) <= capacities[j] * y[j])
        )

    # Create the p-median constraint
    cons_pmedian = m.addCons(quicksum(y[j] for j in range(n_locations)) == n_clusters)

    return m, conss_assignment, conss_capacity, cons_pmedian

def get_simple_instance():
    n = 5
    p = 2
    d = {0: {0: 0, 1: 25, 2: 46, 3: 43, 4: 30}, 1: {1: 0, 2: 22, 3: 20, 4: 22}, 2: {2: 0, 3: 22, 4: 40}, 3: {3: 0, 4: 22}, 4: {4: 0}}
    q = {0: 14, 1: 13, 2: 9, 3: 15, 4: 6}
    Q = {i: 33 for i in range(5)}
    return n, p, d, q, Q

def read_instance_json(path):
    with open(path) as f:
        instance = json.load(f)
    d = {int(k): {int(kk): vv for kk, vv in v.items()} for k, v in instance["d"].items()}
    q = {int(k): v for k, v in instance["q"].items()}
    Q = {int(k): v for k, v in instance["Q"].items()}
    return instance["n"], instance["p"], d, q, Q

if __name__=="__main__":
    n_locations, n_clusters, distances, demands, capacities = read_instance_json("instances/p550-01.json")
    m, *conss = build_model(n_locations, n_clusters, distances, demands, capacities)
    m.optimize()

