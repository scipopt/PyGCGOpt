cdef class VarClassifier:
    """Base class of the Variable Classifier Plugin"""
    cdef public Model model
    cdef public str varclassifiername

    def freeVarClassifier(self):
        '''calls destructor and frees memory of variable classifier'''
        pass

    def classify(self, transformed):
        return {}

cdef SCIP_RETCODE PyVarClassifierFree(SCIP* scip, DEC_VARCLASSIFIER* varclassifier) with gil:
    cdef DEC_CLASSIFIERDATA* varclassifierdata
    varclassifierdata = DECvarClassifierGetData(varclassifier)
    py_varclassifier = <VarClassifier>varclassifierdata
    py_varclassifier.freeVarClassifier()
    Py_DECREF(py_varclassifier)
    return SCIP_OKAY

cdef SCIP_RETCODE PyVarClassifierClassify(SCIP* scip, DEC_VARCLASSIFIER* varclassifier, SCIP_Bool transformed) with gil:
    cdef DEC_CLASSIFIERDATA* varclassifierdata
    varclassifierdata = DECvarClassifierGetData(varclassifier)
    py_varclassifier = <VarClassifier>varclassifierdata
    py_varclassifier.classify(transformed)
    return SCIP_OKAY