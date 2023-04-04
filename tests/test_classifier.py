from pygcgopt import Model, ConsClassifier, VarClassifier, CONS_DECOMPINFO, VAR_DECOMPINFO, quicksum
import pytest

class PyConsClassifier(ConsClassifier):
    def classify(self, conss, partition):
      partition.addClass("demand", "demand conss", CONS_DECOMPINFO.PY_ONLY_PRICING)
      partition.addClass("others", "other conss", CONS_DECOMPINFO.PY_ONLY_MASTER)

      for i in conss:
         if i.name[0:6] == "demand":
            partition.assignConsToClass(i, 0)
         else:
            partition.assignConsToClass(i, 1)

class PyVarClassifier(VarClassifier):
    def classify(self, vars, partition):
      partition.addClass("x", "first x", VAR_DECOMPINFO.PY_BLOCK)
      partition.addClass("others", "first other", VAR_DECOMPINFO.PY_MASTER)

      for i in vars:
         if i.name[0] == "x":
            partition.assignVarToClass(i, 0)
         else:
            partition.assignVarToClass(i, 1)

# taken from https://github.com/scipopt/PySCIPOpt/blob/master/examples/unfinished/cutstock.py
def cuttingStockExample1():
    B = 110                   # roll width (bin size)
    w = [20,45,50,55,75]      # width (size) of orders (items)
    q = [48,35,24,10,8]       # quantitiy of orders
    return w, q, B

def cuttingStockExample2():
    B = 9                     # roll width (bin size)
    w = [2,3,4,5,6,7,8]       # width (size) of orders (items)
    q = [4,2,6,6,2,2,2]       # quantitiy of orders
    return w, q, B

def createCuttingStockModel(width, widths, quantities, ndemands, nrolls):
    model = Model("cutting-stock")

    x = {}
    y = {}

    demands = {}
    capacities = {}

    for i in range(ndemands):
        for j in range(nrolls):
            x[i, j] = model.addVar(vtype="I", name="x_%s_%s"%(i,j))

    for j in range(nrolls):
        y[j] = model.addVar(vtype="B", name="y_%s"%j)

    for i in range(ndemands):
        demands[i] = model.addCons(quicksum(x[i, j] for j in range(nrolls)) >= quantities[i], "demand_%s" % i)
   
    for j in range(nrolls):
        capacities[j] = model.addCons(quicksum(widths[i] * x[i, j] for i in range(ndemands)) <= width * y[j], "capacity_%s" % j)

    model.setObjective(quicksum(y[j] for j in range(nrolls)), "minimize")

    return model, demands, capacities, x, y

@pytest.mark.parametrize("cuttingstockexample", [
    cuttingStockExample1,
    cuttingStockExample2
])
def test_classify(cuttingstockexample):
    widths, quantities, width = cuttingstockexample()
    model, demands, capacities, x, y = createCuttingStockModel(width, widths, quantities, len(widths), sum(quantities))

    for i in model.listConsClassifiers():
        model.setBoolParam("detection/classification/consclassifier/{}/enabled".format(i), False)
    for i in model.listVarClassifiers():
        model.setBoolParam("detection/classification/varclassifier/{}/enabled".format(i), False)

    pyConsClassifier = PyConsClassifier()
    model.includeConsClassifier(pyConsClassifier, "pyconsclassifier", "classify by constraint name")
    assert "pyconsclassifier" in model.listConsClassifiers()

    pyVarClassifier = PyVarClassifier()
    model.includeVarClassifier(pyVarClassifier, "pyvarclassifier", "classify by variable name")
    assert "pyvarclassifier" in model.listVarClassifiers()

    model.setIntParam("presolving/maxrounds", 0)

    model.detect()

    detprobdata = model.getDetprobdataOrig()
    assert detprobdata.getNVarPartitions() == 1
    assert detprobdata.getNConsPartitions() == 1

    # asserts for varpart
    varpartition = detprobdata.getVarPartitions()[0]
    assert varpartition.getName() == "pyvarclassifier"
    assert varpartition.getNClasses() == 2

    # asserts of first class of varpart
    assert varpartition.getClassName(0) == "x"
    assert varpartition.getClassDescription(0) == "first x"
    assert varpartition.getClassDecompInfo(0) == VAR_DECOMPINFO.PY_BLOCK
    assert varpartition.getNVarsOfClasses()[0] == len(widths) * sum(quantities)

    for i in range(len(widths)):
        for j in range(sum(quantities)):
            assert varpartition.getClassOfVar(x[i, j]) == 0
            assert varpartition.getClassNameOfVar(x[i, j]) == "x"

    # asserts of second class of varpart
    assert varpartition.getClassName(1) == "others"
    assert varpartition.getClassDescription(1) == "first other"
    assert varpartition.getClassDecompInfo(1) == VAR_DECOMPINFO.PY_MASTER
    assert varpartition.getNVarsOfClasses()[1] == sum(quantities)

    for j in range(sum(quantities)):
        assert varpartition.getClassOfVar(y[j]) == 1
        assert varpartition.getClassNameOfVar(y[j]) == "others"

    # asserts for conspart
    conspartition = detprobdata.getConsPartitions()[0]
    assert conspartition.getName() == "pyconsclassifier"
    assert conspartition.getNClasses() == 2

    # asserts of first class of conspart
    assert conspartition.getClassName(0) == "demand"
    assert conspartition.getClassDescription(0) == "demand conss"
    assert conspartition.getClassDecompInfo(0) == CONS_DECOMPINFO.PY_ONLY_PRICING
    assert conspartition.getNConssOfClasses()[0] == len(demands)

    for i in range(len(demands)):
        assert conspartition.getClassOfCons(demands[i]) == 0
        assert conspartition.getClassNameOfCons(demands[i]) == "demand"

    # asserts of second class of conspart
    assert conspartition.getClassName(1) == "others"
    assert conspartition.getClassDescription(1) == "other conss"
    assert conspartition.getClassDecompInfo(1) == CONS_DECOMPINFO.PY_ONLY_MASTER
    assert conspartition.getNConssOfClasses()[1] == len(capacities)

    for j in range(len(capacities)):
        assert conspartition.getClassOfCons(capacities[j]) == 1
        assert conspartition.getClassNameOfCons(capacities[j]) == "others"
