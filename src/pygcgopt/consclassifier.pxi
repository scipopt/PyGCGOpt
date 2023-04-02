cdef class ConsClassifier:
    """Base class of the Constraint Classifier Plugin"""
    cdef public Model model
    cdef public str consclassifiername

    def freeConsClassifier(self):
        '''calls destructor and frees memory of constraint classifier'''
        pass

    def classify(self, transformed):
        return {}

cdef SCIP_RETCODE PyConsClassifierFree(SCIP* scip, GCG_CONSCLASSIFIER* consclassifier) with gil:
    cdef GCG_CLASSIFIERDATA* consclassifierdata
    consclassifierdata = GCGconsClassifierGetData(consclassifier)
    py_consclassifier = <ConsClassifier>consclassifierdata
    py_consclassifier.freeConsClassifier()
    Py_DECREF(py_consclassifier)
    return SCIP_OKAY

cdef SCIP_RETCODE PyConsClassifierClassify(SCIP* scip, GCG_CONSCLASSIFIER* consclassifier, SCIP_Bool transformed) with gil:
    cdef GCG_CLASSIFIERDATA* consclassifierdata
    consclassifierdata = GCGconsClassifierGetData(consclassifier)
    py_consclassifier = <ConsClassifier>consclassifierdata
    py_consclassifier.classify(transformed)
    return SCIP_OKAY