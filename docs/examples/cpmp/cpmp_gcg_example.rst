Decomposing the Capacitated p-Median Problem
============================================
The capacitated p-median problem (CPMP) is a well studied problem. Here is a compact formulation of problem:

.. math::

  z_{IP}^{*} = \min \sum_{i=1}^{n} \sum_{j=1}^{n} d_{ij} x_{ij}

There are many different approaches to solve the CPMP. As it has a structure, we can try to solve using a
Branch-Cut-and-Price approach. For that we want to use the Python interface of GCG. We will consider three
use-cases: (1) The Automatic Mode, (2) Exploring different Decompositions, and (3) Building a custom Decomposition.

To follow along with this tutorial interactively, please download the Jupyter notebook and all referenced resources
using the following link: :download:`cpmp_complete.zip<cpmp_gcg_example_files/cpmp_complete.zip>`

Reading in the Instance
_______________________

One of the goals of this tutorial is to demonstrate how a problem modeled in Python can be handed to GCG. Please note that it is also possible to read in a problem in a standard format just like in PySCIPOpt using ``Model.readProb()``.

In this example, we will be using the rather small CPMP instance :download:`p550-01<instances/p550-01.cpmp>`. The file is in a custom file format. We provide the ``reader_cpmp`` module in the resources zip which reads in the instance and returns the input data.

Execute the following block in order to read in the example instance.

.. code-block:: python

  import reader_cpmp

  n_locations, n_clusters, distances, demands, capacities = reader_cpmp.read_instance("instances/p550-01.cpmp")

Setting up the Model
____________________

Now, we want to build the model based on the above formulation. Please note that this part is *not* specific to GCG but is *almost* identical to how one would build the same model with PySCIPOpt. In fact, the only difference is that we import and instanciate ``GCGModel`` instead of ``Model``. In technical terms, ``GCGModel`` is a subclass of ``Model`` and, therefore, you can use all methods of ``Model`` to build your problem.

In order to recreate the model multiple times during this example, we create a method that will return the model. The method also returns the different constraints added to the model grouped by type. This will be important later in use-case 3.

.. code-block:: python

  from pyscipopt import GCGModel, quicksum

  def build_model(n_locations, n_clusters, distances, demands, capacities):
      m = GCGModel()

      m.printVersion()
      m.redirectOutput()

      m.setMinimize()

      x = {}
      y = {}

      for j in range(1, n_locations + 1):
          y[j] = m.addVar(f"y_{j}", vtype="B")
          for i in range(1, n_locations + 1):
              x[i, j] = m.addVar(f"x_{i}_{j}", vtype="B", obj=distances[i][j])

      # Hold different constraint types
      conss_assignment = []
      conss_capacity = []
      cons_pmedian = None

      # Create the assignment constraints
      for i in range(1, n_locations + 1):
          conss_assignment.append(
              m.addCons(quicksum(x[i, j] for j in range(1, n_locations + 1)) == 1)
          )

      # Create the capacity constraints
      for j in range(1, n_locations + 1):
          conss_capacity.append(
              m.addCons(quicksum(demands[i] * x[i, j] for i in range(1, n_locations + 1)) <= capacities[j] * y[j])
          )

      # Create the p-median constraint
      cons_pmedian = m.addCons(quicksum(y[j] for j in range(1, n_locations + 1)) == n_clusters)

      return m, conss_assignment, conss_capacity, cons_pmedian

Use-Case 1: The Automatic Mode
______________________________

With the model built, we can now simply call the ``optimize()`` function and let GCG do its magic.

.. code-block:: python

  m, *conss = build_model(n_locations, n_clusters, distances, demands, capacities)
  m.optimize()

.. code-block::

  presolving:
      (round 1, exhaustive) 0 del vars, 0 del conss, 0 add conss, 0 chg bounds, 0 chg sides, 0 chg coeffs, 100 upgd conss, 0 impls, 50 clqs
         (0.0s) probing: 51/2550 (2.0%) - 0 fixings, 0 aggregations, 0 implications, 0 bound changes
         (0.0s) probing aborted: 50/50 successive totally useless probings
      presolving (2 rounds: 2 fast, 2 medium, 2 exhaustive):
       0 deleted vars, 0 deleted constraints, 0 added constraints, 0 tightened bounds, 0 added holes, 0 changed sides, 0 changed coefficients
       0 implications, 2550 cliques
      presolved problem has 2550 variables (2550 bin, 0 int, 0 impl, 0 cont) and 101 constraints
           50 constraints of type <knapsack>
           50 constraints of type <setppc>
            1 constraints of type <linear>
      transformed objective value is always integral (scale: 1)
      Presolving Time: 0.02
       Consclassifier "nonzeros" yields a classification with 2  different constraint classes
       Consclassifier "constypes" yields a classification with 3 different constraint classes
       Consclassifier "constypes according to miplib" yields a classification with 3 different constraint classes
       Conspartition "constypes according to miplib" is not considered since it offers the same structure as "constypes" conspartition
       Varclassifier "vartypes" yields a classification with 1 different variable classes
       Varclassifier "varobjvals" yields a classification with 116 different variable classes
       Varclassifier "varobjvalsigns" yields a classification with 2 different variable classes
       the current varclass distribution includes 116 classes but only 18 are allowed for GCGconshdlrDecompCalcCandidatesNBlocks()
       in dec_consclass: there are 2 different constraint classes
        the current constraint classifier "nonzeros" consists of 2 different classes
        the current constraint classifier "constypes" consists of 3 different classes
       dec_consclass found 10 new partialdecs
      dec_densemasterconss found 1 new partialdec
      dec_neighborhoodmaster found 1 new partialdec
      POSTPROCESSING of decompositions. Added 0 new decomps.
      Found 11 finished decompositions.
      Measured running time per detector:
      Detector consclass                 worked on        7 finished decompositions and took a total time of      0.000
      Detector neighborhoodmaster        worked on        1 finished decompositions and took a total time of      0.000
      Detector connectedbase             worked on       10 finished decompositions and took a total time of      0.006
      Detector varclass                  worked on        2 finished decompositions and took a total time of      0.001
      Detection Time: 0.01

      A Dantzig-Wolfe reformulation is applied to solve the original problem.
      Chosen structure has 50 blocks and 51 linking constraints.
      This decomposition has a maxwhite score of 0.485149.
      Matrix has 50 blocks, using 50 pricing problems.

        time | node  | left  |SLP iter|MLP iter|LP it/n| mem |mdpt |ovars|mvars|ocons|mcons|mcuts|  dualbound   | primalbound  |  deg   |  gap
      p  0.2s|     1 |     0 |      0 |      0 |     - |  31M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 2.711000e+03 |   --   |    Inf
      p  0.3s|     1 |     0 |      0 |      0 |     - |  30M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 1.785000e+03 |   --   |    Inf
      p  0.3s|     1 |     0 |      0 |      0 |     - |  30M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 1.102000e+03 |   --   |    Inf
         0.3s|     1 |     0 |      0 |      0 |     - |  30M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 1.102000e+03 |   --   |    Inf
         0.3s|     1 |     0 |      0 |      0 |     - |  31M|   0 |2550 |  50 | 102 | 102 |   0 | 0.000000e+00 | 1.102000e+03 |   0.00%|    Inf
         0.3s|     1 |     0 |      0 |      0 |     - |  31M|   0 |2550 | 150 | 102 | 102 |   0 | 0.000000e+00 | 1.102000e+03 |   0.00%|    Inf
         0.5s|     1 |     0 |   2081 |   2081 |     - |  38M|   0 |2550 |1699 | 102 | 102 |   0 | 7.050000e+02 | 1.102000e+03 |  31.56%|  56.31%
      X  1.1s|     1 |     0 |   2747 |   2747 |     - |  50M|   0 |2550 |1749 | 102 | 102 |   0 | 7.050000e+02 | 8.210000e+02 |  31.56%|  16.45%
      Y  1.1s|     1 |     0 |   2747 |   2747 |     - |  50M|   0 |2550 |1799 | 102 | 102 |   0 | 7.050000e+02 | 7.580000e+02 |  31.56%|   7.52%
         1.2s|     1 |     0 |   2747 |   2747 |     - |  50M|   0 |2550 |1799 | 102 | 102 |   0 | 7.050000e+02 | 7.580000e+02 |  31.56%|   7.52%
         1.2s|     1 |     0 |   2747 |   2747 |     - |  50M|   0 |2550 |1799 | 102 | 102 |   0 | 7.050000e+02 | 7.580000e+02 |  31.56%|   7.52%
         1.2s|     1 |     2 |   2747 |   2747 |     - |  50M|   0 |2550 |1799 | 102 | 102 |   0 | 7.050000e+02 | 7.580000e+02 |  31.56%|   7.52%
      *r 1.3s|     6 |     5 |   3484 |   3484 | 147.4 |  51M|   5 |2550 |2016 | 102 | 102 |   0 | 7.050000e+02 | 7.300000e+02 |   --   |   3.55%
      *r 1.3s|     6 |     5 |   3502 |   3502 | 151.0 |  51M|   5 |2550 |2028 | 102 | 102 |   0 | 7.050000e+02 | 7.170000e+02 |   --   |   1.70%
      *r 1.9s|    18 |     6 |   9474 |   9474 | 395.7 |  54M|   5 |2550 |2372 | 102 | 102 |   0 | 7.093000e+02 | 7.130000e+02 |   --   |   0.52%

      SCIP Status        : problem is solved [optimal solution found]
      Solving Time (sec) : 2.79
      Solving Nodes      : 28
      Primal Bound       : +7.13000000000000e+02 (10 solutions)
      Dual Bound         : +7.13000000000000e+02
      Gap                : 0.00 %

Use-Case 2: Exploring different Decompositions
______________________________________________

Above, we have seen GCG in its fully automatic mode. If we want to dig deeper, we can inspect the different decompositions that GCG detects. For that, we recreate the model and manually execute ``presolve()`` and ``detect()`` for the model. At this stage it is possible to list and visualize the found decompositions.

.. code-block:: python

  m, *conss = build_model(n_locations, n_clusters, distances, demands, capacities)
  m.presolve()
  m.detect()

  decomps = m.listDecompositions()

  print("GCG found {} finnished decompositions.".format(len(decomps)))

.. code-block::

  presolving:
      (round 1, exhaustive) 0 del vars, 0 del conss, 0 add conss, 0 chg bounds, 0 chg sides, 0 chg coeffs, 100 upgd conss, 0 impls, 50 clqs
         (0.0s) probing: 51/2550 (2.0%) - 0 fixings, 0 aggregations, 0 implications, 0 bound changes
         (0.0s) probing aborted: 50/50 successive totally useless probings
      presolving (2 rounds: 2 fast, 2 medium, 2 exhaustive):
       0 deleted vars, 0 deleted constraints, 0 added constraints, 0 tightened bounds, 0 added holes, 0 changed sides, 0 changed coefficients
       0 implications, 2550 cliques
      presolved problem has 2550 variables (2550 bin, 0 int, 0 impl, 0 cont) and 101 constraints
           50 constraints of type <knapsack>
           50 constraints of type <setppc>
            1 constraints of type <linear>
      transformed objective value is always integral (scale: 1)
      Presolving Time: 0.02
      starting detection
       Consclassifier "nonzeros" yields a classification with 2  different constraint classes
       Consclassifier "constypes" yields a classification with 3 different constraint classes
       Consclassifier "constypes according to miplib" yields a classification with 3 different constraint classes
       Conspartition "constypes according to miplib" is not considered since it offers the same structure as "constypes" conspartition
       Varclassifier "vartypes" yields a classification with 1 different variable classes
       Varclassifier "varobjvals" yields a classification with 116 different variable classes
       Varclassifier "varobjvalsigns" yields a classification with 2 different variable classes
       the current varclass distribution includes 116 classes but only 18 are allowed for GCGconshdlrDecompCalcCandidatesNBlocks()
       in dec_consclass: there are 2 different constraint classes
        the current constraint classifier "nonzeros" consists of 2 different classes
        the current constraint classifier "constypes" consists of 3 different classes
       dec_consclass found 10 new partialdecs
      dec_densemasterconss found 1 new partialdec
      dec_neighborhoodmaster found 1 new partialdec
      POSTPROCESSING of decompositions. Added 0 new decomps.
      Found 11 finished decompositions.
      Measured running time per detector:
      Detector consclass                 worked on        7 finished decompositions and took a total time of      0.000
      Detector neighborhoodmaster        worked on        1 finished decompositions and took a total time of      0.000
      Detector connectedbase             worked on       10 finished decompositions and took a total time of      0.006
      Detector varclass                  worked on        2 finished decompositions and took a total time of      0.001
      Detection Time: 0.01
      GCG found 11 finnished decompositions.

Inspecting Decompositions
^^^^^^^^^^^^^^^^^^^^^^^^^

The call to ``listDecompositions()`` returns a list of ``PartialDecomposition`` objects. We can print a decomposition using the Python ``print()`` function to get a summary or access different fields directly.

For a full overview of available methods, take a look at the online documentation for the ``PartialDecomposition`` class, or execute ``help(d)`` where ``d`` is a decomposition object.

.. code-block:: python

  print(decomps)

  d = decomps[2]

  print(
      f"Decomp scores: {d.classicScore=:.04f}, {d.borderAreaScore=:.04f}, {d.maxWhiteScore=:.04f}, {d.maxForWhiteScore=:.04f}"
  )

.. code-block::

  [<PartialDecomposition: nBlocks=0, nMasterConss=101, nMasterVars=2550, nLinkingVars=0, maxForWhiteScore=0.0>, <PartialDecomposition: nBlocks=1, nMasterConss=0, nMasterVars=0, nLinkingVars=0, maxForWhiteScore=0.0>, <PartialDecomposition: nBlocks=50, nMasterConss=51, nMasterVars=0, nLinkingVars=0, maxForWhiteScore=0.48514851485148514>, <PartialDecomposition: nBlocks=51, nMasterConss=50, nMasterVars=0, nLinkingVars=0, maxForWhiteScore=0.49504950495049505>, <PartialDecomposition: nBlocks=1, nMasterConss=50, nMasterVars=0, nLinkingVars=0, maxForWhiteScore=0.0>, <PartialDecomposition: nBlocks=1, nMasterConss=100, nMasterVars=2500, nLinkingVars=0, maxForWhiteScore=0.009706853038245034>, <PartialDecomposition: nBlocks=1, nMasterConss=1, nMasterVars=0, nLinkingVars=0, maxForWhiteScore=0.0>, <PartialDecomposition: nBlocks=50, nMasterConss=51, nMasterVars=50, nLinkingVars=0, maxForWhiteScore=0.4853426519122501>, <PartialDecomposition: nBlocks=1, nMasterConss=1, nMasterVars=0, nLinkingVars=0, maxForWhiteScore=0.0>, <PartialDecomposition: nBlocks=101, nMasterConss=0, nMasterVars=0, nLinkingVars=2550, maxForWhiteScore=0.019291161956034086>, <PartialDecomposition: nBlocks=0, nMasterConss=101, nMasterVars=100, nLinkingVars=2450, maxForWhiteScore=0.0>]
  Decomp scores: d.classicScore=-1.0000, d.borderAreaScore=-1.0000, d.maxWhiteScore=0.4851, d.maxForWhiteScore=0.4851

Visualizing Decompositions
^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition, GCG can create graphical visualizations of decompositions. They can easily be displayed in a Jupyter nodebook like so:

.. code-block:: python

  d

![svg](cpmp_gcg_example_files/cpmp_gcg_example_12_0.svg)

Selecting Decompositions
^^^^^^^^^^^^^^^^^^^^^^^^

After this investigation, we decide that we like this particular decomposition. The following code orders GCG to select our decomposition instead of an automatic one. Then, the optimization process is started.

.. code-block:: python

  d.isSelected = True

  m.optimize()

.. code-block::

  A Dantzig-Wolfe reformulation is applied to solve the original problem.
    Chosen structure has 50 blocks and 51 linking constraints.
    This decomposition has a maxwhite score of 0.485149.
    Matrix has 50 blocks, using 50 pricing problems.

      time | node  | left  |SLP iter|MLP iter|LP it/n| mem |mdpt |ovars|mvars|ocons|mcons|mcuts|  dualbound   | primalbound  |  deg   |  gap
    p  0.2s|     1 |     0 |      0 |      0 |     - |  31M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 2.711000e+03 |   --   |    Inf
    p  0.3s|     1 |     0 |      0 |      0 |     - |  30M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 1.785000e+03 |   --   |    Inf
    p  0.3s|     1 |     0 |      0 |      0 |     - |  30M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 1.102000e+03 |   --   |    Inf
       0.3s|     1 |     0 |      0 |      0 |     - |  30M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 1.102000e+03 |   --   |    Inf
       0.3s|     1 |     0 |      0 |      0 |     - |  31M|   0 |2550 |  50 | 102 | 102 |   0 | 0.000000e+00 | 1.102000e+03 |   0.00%|    Inf
       0.3s|     1 |     0 |      0 |      0 |     - |  31M|   0 |2550 | 150 | 102 | 102 |   0 | 0.000000e+00 | 1.102000e+03 |   0.00%|    Inf
       0.5s|     1 |     0 |   2081 |   2081 |     - |  38M|   0 |2550 |1699 | 102 | 102 |   0 | 7.050000e+02 | 1.102000e+03 |  31.56%|  56.31%
    X  1.1s|     1 |     0 |   2747 |   2747 |     - |  50M|   0 |2550 |1749 | 102 | 102 |   0 | 7.050000e+02 | 8.210000e+02 |  31.56%|  16.45%
    Y  1.1s|     1 |     0 |   2747 |   2747 |     - |  50M|   0 |2550 |1799 | 102 | 102 |   0 | 7.050000e+02 | 7.580000e+02 |  31.56%|   7.52%
       1.1s|     1 |     0 |   2747 |   2747 |     - |  50M|   0 |2550 |1799 | 102 | 102 |   0 | 7.050000e+02 | 7.580000e+02 |  31.56%|   7.52%
       1.2s|     1 |     0 |   2747 |   2747 |     - |  50M|   0 |2550 |1799 | 102 | 102 |   0 | 7.050000e+02 | 7.580000e+02 |  31.56%|   7.52%
       1.2s|     1 |     2 |   2747 |   2747 |     - |  50M|   0 |2550 |1799 | 102 | 102 |   0 | 7.050000e+02 | 7.580000e+02 |  31.56%|   7.52%
    *r 1.3s|     6 |     5 |   3484 |   3484 | 147.4 |  51M|   5 |2550 |2016 | 102 | 102 |   0 | 7.050000e+02 | 7.300000e+02 |   --   |   3.55%
    *r 1.3s|     6 |     5 |   3502 |   3502 | 151.0 |  51M|   5 |2550 |2028 | 102 | 102 |   0 | 7.050000e+02 | 7.170000e+02 |   --   |   1.70%
    *r 1.9s|    18 |     6 |   9474 |   9474 | 395.7 |  54M|   5 |2550 |2372 | 102 | 102 |   0 | 7.093000e+02 | 7.130000e+02 |   --   |   0.52%

    SCIP Status        : problem is solved [optimal solution found]
    Solving Time (sec) : 2.76
    Solving Nodes      : 28
    Primal Bound       : +7.13000000000000e+02 (10 solutions)
    Dual Bound         : +7.13000000000000e+02
    Gap                : 0.00 %

Use-Case 3: Building a custom Decomposition
___________________________________________

In the previous use-cases we run GCG with automatically detected decompositions. This is useful if we do not know the exact structure of a model but still want to exploit GCG's decomposition approach.

However, for our model we *do* know its structure. If you take another look at our ``build_model()`` method, you will notice that we store the created constraints in different variables based on their type. This is a typical approach when we want to specify a custom decomposition after building the model using the Python interface.

In the following code, we recreate our model and use the different constraint types fo select constraints for reformulation and constraints for the Dantzig-Wolfe master problem.

.. code-block:: python

  m, conss_assignment, conss_capacity, cons_pmedian = build_model(
      n_locations, n_clusters, distances, demands, capacities
  )

  conss_master = conss_assignment + [cons_pmedian]
  conss_reform = conss_capacity

  pd = m.createPartialDecomposition()
  pd.fixConssToMaster(conss_master)

  for block, c in enumerate(conss_reform):
      pd.fixConsToBlock(c, block)

  m.addPreexistingPartialDecomposition(pd)

  m.optimize()

.. code-block::

  added complete decomp for original problem with 50 blocks and 51 masterconss, 0 linkingvars, 0 mastervars, and max white score of   0.485149
    there is an original decomposition and problem is not presolved yet -> disable presolving and start optimizing (rerun with presolve command before detect command for detecting in presolved problem)
    presolving:
    presolving (0 rounds: 0 fast, 0 medium, 0 exhaustive):
     0 deleted vars, 0 deleted constraints, 0 added constraints, 0 tightened bounds, 0 added holes, 0 changed sides, 0 changed coefficients
     0 implications, 0 cliques
    presolved problem has 2550 variables (2550 bin, 0 int, 0 impl, 0 cont) and 101 constraints
        101 constraints of type <linear>
    transformed objective value is always integral (scale: 1)
    Presolving Time: 0.00
     calculated translation; number of missing constraints: 0; number of other partialdecs: 1
    Preexisting decomposition found. Solution process started...

    A Dantzig-Wolfe reformulation is applied to solve the original problem.
    Chosen structure has 50 blocks and 51 linking constraints.
    This decomposition has a maxwhite score of 0.485149.
    Matrix has 50 blocks, using 50 pricing problems.

      time | node  | left  |SLP iter|MLP iter|LP it/n| mem |mdpt |ovars|mvars|ocons|mcons|mcuts|  dualbound   | primalbound  |  deg   |  gap
    p  0.2s|     1 |     0 |      0 |      0 |     - |  27M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 2.848000e+03 |   --   |    Inf
       0.2s|     1 |     0 |      0 |      0 |     - |  27M|   0 |2550 |   0 | 101 |   0 |   0 | 0.000000e+00 | 2.848000e+03 |   --   |    Inf
       0.2s|     1 |     0 |      0 |      0 |     - |  28M|   0 |2550 |  50 | 102 | 102 |   0 | 0.000000e+00 | 2.848000e+03 |   0.00%|    Inf
       0.2s|     1 |     0 |      0 |      0 |     - |  28M|   0 |2550 |  50 | 102 | 102 |   0 | 0.000000e+00 | 2.848000e+03 |   0.00%|    Inf
       0.4s|     1 |     0 |   1839 |   1839 |     - |  38M|   0 |2550 |2302 | 102 | 102 |   0 | 0.000000e+00 | 2.848000e+03 |  34.27%|    Inf
       0.6s|     1 |     0 |   3939 |   3939 |     - |  40M|   0 |2550 |2933 | 102 | 102 |   0 | 7.024809e+02 | 2.848000e+03 |  27.83%| 305.42%
       0.6s|     1 |     0 |   3939 |   3939 |     - |  40M|   0 |2550 |2933 | 102 | 102 |   0 | 7.050000e+02 | 2.848000e+03 |  27.98%| 303.97%
    X  1.1s|     1 |     0 |   4628 |   4628 |     - |  46M|   0 |2550 |2983 | 102 | 102 |   0 | 7.050000e+02 | 7.950000e+02 |  27.98%|  12.77%
       1.2s|     1 |     0 |   4628 |   4628 |     - |  46M|   0 |2550 |3033 | 102 | 102 |   0 | 7.050000e+02 | 7.950000e+02 |  27.98%|  12.77%
       1.2s|     1 |     0 |   4628 |   4628 |     - |  46M|   0 |2550 |3033 | 102 | 102 |   0 | 7.050000e+02 | 7.950000e+02 |  27.98%|  12.77%
       1.2s|     1 |     2 |   4628 |   4628 |     - |  46M|   0 |2550 |3033 | 102 | 102 |   0 | 7.050000e+02 | 7.950000e+02 |  27.98%|  12.77%
    *r 1.3s|     4 |     3 |   5331 |   5331 | 234.3 |  47M|   3 |2550 |3106 | 102 | 102 |   0 | 7.050000e+02 | 7.140000e+02 |   --   |   1.28%
    *r 1.6s|     9 |     6 |   6230 |   6230 | 200.2 |  48M|   3 |2550 |3310 | 102 | 102 |   0 | 7.097500e+02 | 7.130000e+02 |   --   |   0.46%

    SCIP Status        : problem is solved [optimal solution found]
    Solving Time (sec) : 2.05
    Solving Nodes      : 19
    Primal Bound       : +7.13000000000000e+02 (9 solutions)
    Dual Bound         : +7.13000000000000e+02
    Gap                : 0.00 %

Summary
_______

With that, we have seen the most important features to use GCG as a solver through the Python interface. In a different example, we will take a look at how GCG can be extended using Python code.
