from pygcgopt import GCGModel, quicksum as qs

n_locs = 5
n_clusters = 2
distances = {0: {0: 0, 1: 6, 2: 54, 3: 52, 4: 19}, 1: {0: 6, 1: 0, 2: 28, 3: 75, 4: 61}, 2: {0: 54, 1: 28, 2: 0, 3: 91, 4: 40}, 3: {0: 52, 1: 75, 2: 91, 3: 0, 4: 28}, 4: {0: 19, 1: 61, 2: 40, 3: 28, 4: 0},}
demands = {0: 14, 1: 13, 2: 9, 3: 15, 4: 6}
capacities = {0: 39, 1: 39, 2: 39, 3: 39, 4: 39}

m = GCGModel()
x = {(i, j): m.addVar(f"x_{i}_{j}", vtype="B", obj=distances[i][j]) for i in range(n_locs) for j in range(n_locs)}
y = {j: m.addVar(f"y_{j}", vtype="B") for j in range(n_locs)}

conss_assignment = m.addConss(
  [qs(x[i, j] for j in range(n_locs)) == 1 for i in range(n_locs)])
conss_capacity = m.addConss(
  [qs(demands[i] * x[i, j] for i in range(n_locs)) <= capacities[j] * y[j] for j in range(n_locs)])
cons_pmedian = m.addCons(qs(y[j] for j in range(n_locs)) == n_clusters)

master_conss = conss_assignment + [cons_pmedian]
block_conss = [[cons] for cons in conss_capacity]
m.addDecompositionFromConss(master_conss, *block_conss)

m.optimize()
