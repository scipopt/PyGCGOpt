cdef class VarClassifier:
    """Base class of the Variable Classifier Plugin"""
    cdef public Model model
    cdef public str varclassifiername

    def freeVarClassifier(self):
        '''calls destructor and frees memory of variable classifier'''
        pass

    def classify(self, detprobdata):
        return {}

cdef SCIP_RETCODE PyVarClassifierFree(SCIP* scip, GCG_VARCLASSIFIER* varclassifier) with gil:
    cdef GCG_CLASSIFIERDATA* varclassifierdata
    varclassifierdata = GCGvarClassifierGetData(varclassifier)
    py_varclassifier = <VarClassifier>varclassifierdata
    py_varclassifier.freeVarClassifier()
    Py_DECREF(py_varclassifier)
    return SCIP_OKAY

cdef SCIP_RETCODE PyVarClassifierClassify(SCIP* scip, GCG_VARCLASSIFIER* varclassifier, SCIP_Bool transformed) with gil:
    cdef GCG_CLASSIFIERDATA* varclassifierdata
    varclassifierdata = GCGvarClassifierGetData(varclassifier)
    py_varclassifier = <VarClassifier>varclassifierdata
    if transformed:
        detprobdata = py_varclassifier.model.getDetprobdataPresolved()
    else:
        detprobdata = py_varclassifier.model.getDetprobdataOrig()
    py_varclassifier.classify(detprobdata)
    return SCIP_OKAY
