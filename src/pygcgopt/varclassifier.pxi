cdef class VarClassifier:
    """Base class of the Variable Classifier Plugin"""
    cdef public Model model
    cdef public str name

    def freeVarClassifier(self):
        '''calls destructor and frees memory of variable classifier'''
        pass

    def classify(self, vars, partition):
        return {}

cdef SCIP_RETCODE PyVarClassifierFree(GCG* gcg, GCG_CLSVAR* varclassifier) with gil:
    cdef GCG_CLSVARDATA* varclassifierdata
    varclassifierdata = GCGclsvarGetData(varclassifier)
    py_varclassifier = <VarClassifier>varclassifierdata
    py_varclassifier.freeVarClassifier()
    Py_DECREF(py_varclassifier)
    return SCIP_OKAY

cdef SCIP_RETCODE PyVarClassifierClassify(GCG* gcg, GCG_CLSVAR* varclassifier, SCIP_Bool transformed) with gil:
    cdef GCG_CLSVARDATA* varclassifierdata
    varclassifierdata = GCGclsvarGetData(varclassifier)
    py_varclassifier = <VarClassifier>varclassifierdata
    if transformed:
        detprobdata = py_varclassifier.model.getDetprobdataPresolved()
    else:
        detprobdata = py_varclassifier.model.getDetprobdataOrig()
    vars = detprobdata.getRelevantVars()
    partition = detprobdata.createVarPart(py_varclassifier.name, 0, len(vars))
    py_varclassifier.classify(vars, partition)
    print("Varclassifier {0} yields a classification with {1} different variable classes".format(partition.getName(), partition.getNClasses()))
    detprobdata.addVarPartition(partition)
    return SCIP_OKAY
