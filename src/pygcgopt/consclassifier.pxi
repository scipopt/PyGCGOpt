cdef class ConsClassifier:
    """Base class of the Constraint Classifier Plugin"""
    cdef public Model model
    cdef public str consclassifiername

    def freeConsClassifier(self):
        '''calls destructor and frees memory of constraint classifier'''
        pass

    def classify(self, transformed):
        return {}

cdef SCIP_RETCODE PyConsClassifierFree(SCIP* scip, DEC_CONSCLASSIFIER* consclassifier) with gil:
    cdef DEC_CLASSIFIERDATA* consclassifierdata
    consclassifierdata = DECconsClassifierGetData(consclassifier)
    py_consclassifier = <ConsClassifier>consclassifierdata
    py_consclassifier.freeConsClassifier()
    Py_DECREF(py_consclassifier)
    return SCIP_OKAY

cdef SCIP_RETCODE PyConsClassifierClassify(SCIP* scip, DEC_CONSCLASSIFIER* consclassifier, SCIP_Bool transformed) with gil:
    cdef DEC_CLASSIFIERDATA* consclassifierdata
    consclassifierdata = DECconsClassifierGetData(consclassifier)
    py_consclassifier = <ConsClassifier>consclassifierdata
    py_consclassifier.classify(transformed)
    return SCIP_OKAY