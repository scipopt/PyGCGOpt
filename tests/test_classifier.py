from pygcgopt import Model, VarClassifier, VAR_DECOMPINFO, quicksum

#import pytest

class PyVarClassifier(VarClassifier):
    def classify(self, vars, partition):
      partition.addClass("x", "first x", VAR_DECOMPINFO.PY_BLOCK)
      partition.addClass("other", "first other", VAR_DECOMPINFO.PY_MASTER)

      for i in vars:
         if i.name[0] == "x":
            partition.assignVarToClass(i, 0)
         else:
            partition.assignVarToClass(i, 1)


def CuttingStockExample1():
    """CuttingStockExample1: create toy instance for the cutting stock problem."""
    B = 110            # roll width (bin size)
    w = [20,45,50,55,75]  # width (size) of orders (items)
    q = [48,35,24,10,8]  # quantitiy of orders
    return w, q, B

def CuttingStockExample2():
    """CuttingStockExample2: create toy instance for the cutting stock problem."""
    B = 9            # roll width (bin size)
    w = [2,3,4,5,6,7,8]   # width (size) of orders (items)
    q = [4,2,6,6,2,2,2]  # quantitiy of orders
    return w, q, B

def createCuttingStockModel(width, widths, quantities, ndemands, nrolls):
    model = Model("cutting-stock")

    x = {}
    y = {}

    for i in range(ndemands):
        for j in range(nrolls):
            x[i,j] = model.addVar(vtype="I", name="x_%s_%s"%(i,j))

    for j in range(nrolls):
        y[j] = model.addVar(vtype="B", name="y_%s"%j)

    for i in range(ndemands):
        model.addCons(quicksum(x[i, j] for j in range(nrolls)) >= quantities[i], "Demand_%s" % i)
   
    for j in range(nrolls):
        model.addCons(quicksum(widths[i] * x[i, j] for i in range(ndemands)) <= width * y[j], "Capacity_%s" % j)

    model.setObjective(quicksum(y[j] for j in range(nrolls)), "minimize")

    return model

def test_classify():
    widths, quantities, width = CuttingStockExample1()
    model = createCuttingStockModel(width, widths, quantities, len(widths), sum(quantities))

    for i in model.listVarClassifiers():
        model.setBoolParam("detection/classification/varclassifier/{}/enabled".format(i), False)
    for i in model.listConsClassifiers():
        model.setBoolParam("detection/classification/consclassifier/{}/enabled".format(i), False)
    pyVarClassifier = PyVarClassifier()
    model.includeVarClassifier(pyVarClassifier, "pyvarclassifier", "classify by variable name")
    assert "py-varname" in model.listVarClassifiers()

    model.setIntParam("presolving/maxrounds", 0)

    model.detect()

    detprobdata = model.getDetprobdataOrig()
    assert detprobdata.getNVarPartitions() == 1

if __name__ == "__main__":
    test_score()
