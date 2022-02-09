cdef class ConsClassifier:
    """Base class of the Constraint Classifier Plugin"""
    cdef public SCIPModel model
    cdef public str consclassifiername

    def freeConsClassifier(self):
        '''calls destructor and frees memory of constraint classifier'''
        pass

    def classify(self, transformed):
        return {}

cdef ConsClassifier get_py_consclassifier(DEC_CONSCLASSIFIER* consclassifier):
    cdef DEC_CLASSIFIERDATA* consclassifierdata
    consclassifierdata = DECconsClassifierGetData(consclassifier)
    py_consclassifier = <ConsClassifier>consclassifierdata
    return py_consclassifier

cdef SCIP_RETCODE PyConsClassifierFree(SCIP* scip, DEC_CONSCLASSIFIER* consclassifier):
    py_consclassifier = get_py_consclassifier(consclassifier)
    py_consclassifier.freeConsClassifier()
    Py_DECREF(py_consclassifier)
    return SCIP_OKAY

cdef SCIP_RETCODE PyConsClassifierClassify(SCIP* scip, DEC_CONSCLASSIFIER* consclassifier, SCIP_Bool transformed):
    py_consclassifier = get_py_consclassifier(consclassifier)
    py_consclassifier.classify(transformed)
    return SCIP_OKAY