Exploring DW Dual Bounds
========================

A main feature of GCG is its plugin interface which can be used to plug into different parts of the solving process. In this tutorial we will take a look at how to implement a simple detector plugin. In GCG, detectors are executed before the solving process to find structures in the model passed to the solver. For more details in detectors, please consult GCG's documentation.

The example used in this tutorial is based on a paper by Witt et al. We will built a detector in Python that detects all possible decompositions for a given model. After that, we will run a simple experiment in collect statistics using the Python interface. To follow along with this tutorial interactively, please download the Jupyter notebook and all referenced resources
using the following link: :download:`alldecomps_complete.zip<../cpmp/cpmp_gcg_example_files/cpmp_complete.zip>`

Preparation
___________

We start with adding all neccessary imports as well as a utility function to return the powerset of a list. This will be usefull later since reformulating each set of constraints in the powerset is equivalent to obtaining all decompositions.

.. code-block:: python

  from pygcgopt import GCGModel, Detector, SCIP_PARAMSETTING, SCIP_STAGE

  from itertools import chain, combinations
  import json
  from datetime import datetime
  from pathlib import Path
  import sys

  def powerset(iterable):
      """Returns the powerset of the passed iterable."""
      s = list(iterable)
      return map(set, chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))

  list(powerset([0,1,2]))

::

  [set(), {0}, {1}, {2}, {0, 1}, {0, 2}, {1, 2}, {0, 1, 2}]

The Detector
____________

To implement a detector for GCG in Python, you need to subclass the ``Detector`` class. It has several methods that can be overriden to plugin at different stages of the detector loop. For more details on the other methods, please consult the documentation of the Python interface. For general information on the detector loop, please consult GCG's documentation. The callbacks of the Python interface match those of the C interface.

Without further ado, here is the code for the detector. We will take a look at the different parts of it afterwards.

.. code-block:: python

  class AllDecompsDetector(Detector):
      def __init__(self, n_conss):
          super().__init__()

          self.conss_powerset = powerset(range(n_conss))
          self.iteration_idx = 0

      def propagatePartialdec(self, detprobdata, workonpartialdec):
          all_conss = sorted(workonpartialdec.getOpenconss(), key=lambda cons: cons.name)

          new_decs = []
          for reform_conss_idx in self.conss_powerset:
              reform_conss = set(all_conss[i] for i in reform_conss_idx)
              master_conss = set(all_conss) - reform_conss
              assert len(master_conss) + len(reform_conss) == len(all_conss)

              new_dec = workonpartialdec.copy()
              new_dec.fixConssToMaster(master_conss)
              new_dec.fixConssToBlock(reform_conss, 0)

              new_decs.append(new_dec)

              # print(f"Produced decomposition at index {self.iteration_idx}: {reform_conss}")
              self.last_decomp = list(c.name for c in reform_conss)

              self.iteration_idx += 1

              # The unconditional break is intentional, read on to find out why.
              break


          return {
              "newpartialdecs": [new_dec]
          }

The above class is already enough to define our detector!

Let's walk it through:

In the constructor, we compute the powerset based on the number of constraints our model has. Each element of the ``conss_powerset`` list will be a *set* of constraints to reformulate.

The ``propagatePartialdec`` method is one of three methods that can be overriden for a detector. It receives two arguments ``detprobdata`` of type ``DetProbData`` and ``workonpartialdec`` of type ``PartialDecomposition``. The first argument contains generic information accumulated during the detection and classification process. The second argument contains a partial decomposition that our detector *may* use to derive new (partial) decompositions. In our case, the passed partial decomposition will always be empty due to how we will setup the experiment later.

Our detector gets all constraints from the partial decomposition and then iterates over the previously generated powerset of constraints. We select the constraints to reformulate and the constraints to remain in the master problem. Then, we copy the partial decomposition and assign the master and reformulation constraints. In the end, we return all new partial decompositions.

Note: You may wonder why we use a for-loop to iterate over the powerset but exit after one iteration. Down below, we will create a fresh ``Model`` for *every* decomposition and assign the *same* instance of the detector. GCG will call the ``propagatePartialdec`` method once for each model and our code will return a different decomposition each time because ``conss_powerset`` is a Python iterable object. Thus, the iteration will *not* start over but ascend every time.

Experiment
__________

Now that we have our detector, we want to use it to replicate the study referanced above.

Setting up the Model
^^^^^^^^^^^^^^^^^^^^

The goal is to solve the LP relaxation of the RMP of *every* reformulation. For that, we need to disable cutting planes, limit to number of nodes, and prevent GCG from aborting pricing *before* the RMP is optimal. We create a function to create a fresh model and make these settings.

.. note::
  The parameter ``limits/nodes`` is not a GCG parameter but a SCIP parameter. We can set using the appropriate method of PySCIPOpt's ``Model`` class. You can set any SCIP or GCG parameter in a similar manner.

.. code-block:: python

  def init_model():
      m = GCGModel()

      for det in m.listDetectors():
          m.setDetectorEnabled(det, False)
          m.setDetectorFinishingEnabled(det, False)
          m.setDetectorPostprocessingEnabled(det, False)

      m.setGCGSeparating(SCIP_PARAMSETTING.OFF)
      m.setLongintParam("limits/nodes", 1)
      m.setBoolParam("pricing/masterpricer/abortpricingint", False)

      return m

Running the Experiment
^^^^^^^^^^^^^^^^^^^^^^

Now comes the fun part, we can run our little experiment!

While the following code might look terrifying at first, you will notice that most of it is just boilerplate code to setup logging and to collect statistics. The important line is the one containing the call to ``m.includeDetector()``. This registers our detector to GCG which will use it in addition to its predefined detectors. In our case, we disabled the other detectors when creating the model and, therefore, GCG will only use our detector.

The general procedure is as follows:
  1. Create a fresh model using ``init_model()``
  2. Include our detector
  3. Read in our problem instance
  4. Detect and solve the problem (we explicitely call ``detect()`` to avoid presolving the problem, see documentation for details)
  5. Collect and store statistics

Running the experiment will take a few minutes, don't panic!

.. code-block:: python

  results_dir = Path("results/").joinpath(problem_name)
  results_dir.mkdir(parents=True, exist_ok=True)
  logs_dir = results_dir.joinpath("log")
  logs_dir.mkdir(exist_ok=True)

  current_timestamp = datetime.utcnow().strftime("%Y-%m-%dT%H%M%S")
  results_file = results_dir.joinpath(f"result_{problem_name}_{current_timestamp}.jsonl").resolve()

  for i in range(2**n_conss):
      log_path = logs_dir.joinpath(f"log_{problem_name}_{current_timestamp}_idx_{i:06}.log").resolve().as_posix()

      m = init_model()
      m.setLogfile(log_path)
      m.includeDetector(all_decomps_detector, "all_decomps_detector", "a", "Detects the power set of constraints")
      m.readProblem(problem_path)

      m.detect()

      assert(len(m.listDecompositions()) == 1)

      m.optimize()

      mp = m.getMasterProb()

      result = {
          "reformulation_constraints": all_decomps_detector.last_decomp,
          "iteration_idx": i,
          "dual_bound": m.getDualbound(),
          "total_time": m.getTotalTime(),
          "solving_time": m.getSolvingTime(),
          "reading_time": m.getReadingTime(),
          "presolving_time": m.getPresolvingTime(),
          "status": m.getStatus(),
          "log_filename": log_path,
      }

      with results_file.open("a") as f:
          json.dump(result, f)
          f.write('\n')

      m.freeProb()

  print("Finished experiment!")

::

  Finished experiment!

Summary
_______

With that we finished our little experiment and obtained the dual bounds of all possible reformulations. In a next step, we would evaluate the data that is stored in the results file.