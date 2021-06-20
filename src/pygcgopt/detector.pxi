##@file detector.pxi
#@brief Base class of the Detector Plugin
cdef class Detector:
    cdef public Model model
    cdef public str detectorname

    def freeDetector(self):
        '''calls destructor and frees memory of detector'''
        pass

    def initDetector(self):
        '''initializes detector'''
        pass

    def exitDetector(self):
        '''calls exit method of detector'''
        pass

    def propagatePartialdec(self, detprobdata, workonpartialdec):
        '''refines a partial decomposition inside detection loop'''
        return {}

    def finishPartialdec(self, detprobdata, workonpartialdec):
        '''completes a partial decomposition when called in detection loop'''
        return {}

    def postprocessPartialdec(self, detprobdata, workonpartialdec):
        '''postprocess a complete decomposition, called after detection loop'''
        return {}

    def setParamAggressive(self):
        '''called if the detection emphasis setting aggressive is chosen'''
        pass

    def setParamDefault(self):
        '''called if the detection emphasis setting default is chosen'''
        pass

    def setParamFast(self):
        '''called if the detection emphasis setting fast is chosen'''
        pass


cdef Detector get_py_detector(DEC_DETECTOR* detector):
    cdef DEC_DETECTORDATA* detectordata
    detectordata = DECdetectorGetData(detector)
    py_detector = <Detector>detectordata
    return py_detector


cdef tuple get_py_partialdec_detection_data(PARTIALDEC_DETECTION_DATA* partialdecdetectiondata):
    detprobdata = DetProbData.create(partialdecdetectiondata.detprobdata)
    workonpartialdec = PartialDecomposition.create(partialdecdetectiondata.workonpartialdec)
    return detprobdata, workonpartialdec


cdef wrap_detector_callback_result(Detector detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result, object result_dict, double detection_time):
    if result_dict is None or type(result_dict) is not dict:
        raise TypeError("Detector callback returned '{}' of type '{}'. Expected a dictionary with optional keys 'newpartialdecs' and 'result'.".format(result_dict, type(result_dict)))

    py_newpartialdecs = result_dict.get("newpartialdecs", [])
    if py_newpartialdecs is None or type(py_newpartialdecs) is not list:
        raise TypeError("The value of the key 'newpartialdecs' is '{}' of type '{}'. Expected a list of 'PartialDecomposition' objects.".format(py_newpartialdecs, type(py_newpartialdecs)))

    nnewpartialdecs = len(py_newpartialdecs)
    cdef PARTIALDECOMP** newpartialdecs = <PARTIALDECOMP**>malloc(nnewpartialdecs * sizeof(PARTIALDECOMP*))
    for i in range(nnewpartialdecs):
        py_newpartialdecs[i].sort()
        # TODO: It would be nice if the user could set a custom chain info
        py_newpartialdecs[i].addDetectorChainInfo(detector.detectorname)
        py_newpartialdecs[i].addClockTime(detection_time / nnewpartialdecs)
        newpartialdecs[i] = (<PartialDecomposition>py_newpartialdecs[i]).thisptr
    partialdecdetectiondata.newpartialdecs = newpartialdecs
    partialdecdetectiondata.nnewpartialdecs = nnewpartialdecs
    partialdecdetectiondata.detectiontime = detection_time
    result[0] = result_dict.get("result", <SCIP_RESULT>result[0])


cdef SCIP_RETCODE PyDetectorFree (SCIP* scip, DEC_DETECTOR* detector):
    py_detector = get_py_detector(detector)
    py_detector.freeDetector()
    Py_DECREF(py_detector)
    return SCIP_OKAY

cdef SCIP_RETCODE PyDetectorInit (SCIP* scip, DEC_DETECTOR* detector):
    py_detector = get_py_detector(detector)
    py_detector.initDetector()
    return SCIP_OKAY

cdef SCIP_RETCODE PyDetectorExit (SCIP* scip, DEC_DETECTOR* detector):
    py_detector = get_py_detector(detector)
    py_detector.exitDetector()
    return SCIP_OKAY

cdef SCIP_RETCODE PyDetectorPropagatePartialdec (SCIP* scip, DEC_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result):
    cdef SCIP_CLOCK* clock = start_new_clock(scip)

    py_detector = get_py_detector(detector)
    detprobdata, workonpartialdec = get_py_partialdec_detection_data(partialdecdetectiondata)

    result_dict = py_detector.propagatePartialdec(detprobdata, workonpartialdec)

    cdef double detection_time = stop_and_free_clock(scip, clock)
    wrap_detector_callback_result(py_detector, partialdecdetectiondata, result, result_dict, detection_time)

    return SCIP_OKAY

cdef SCIP_RETCODE PyDetectorFinishPartialdec (SCIP* scip, DEC_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result):
    cdef SCIP_CLOCK* clock = start_new_clock(scip)

    py_detector = get_py_detector(detector)
    detprobdata, workonpartialdec = get_py_partialdec_detection_data(partialdecdetectiondata)

    result_dict = py_detector.finishPartialdec(detprobdata, workonpartialdec)

    cdef double detection_time = stop_and_free_clock(scip, clock)
    wrap_detector_callback_result(py_detector, partialdecdetectiondata, result, result_dict, detection_time)

    return SCIP_OKAY

cdef SCIP_RETCODE PyDetectorPostprocessPartialdec (SCIP* scip, DEC_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result):
    cdef SCIP_CLOCK* clock = start_new_clock(scip)

    py_detector = get_py_detector(detector)
    detprobdata, workonpartialdec = get_py_partialdec_detection_data(partialdecdetectiondata)

    result_dict = py_detector.postprocessPartialdec(detprobdata, workonpartialdec)

    cdef double detection_time = stop_and_free_clock(scip, clock)
    wrap_detector_callback_result(py_detector, partialdecdetectiondata, result, result_dict, detection_time)

    return SCIP_OKAY

cdef SCIP_RETCODE PyDetectorSetParamAggressive (SCIP* scip, DEC_DETECTOR* detector, SCIP_RESULT* result):
    py_detector = get_py_detector(detector)
    result_dict = py_detector.setParamAggressive() or {}
    result[0] = result_dict.get("result", <SCIP_RESULT>result[0])
    return SCIP_OKAY

cdef SCIP_RETCODE PyDetectorSetParamDefault (SCIP* scip, DEC_DETECTOR* detector, SCIP_RESULT* result):
    py_detector = get_py_detector(detector)
    result_dict = py_detector.setParamDefault() or {}
    result[0] = result_dict.get("result", <SCIP_RESULT>result[0])
    return SCIP_OKAY

cdef SCIP_RETCODE PyDetectorSetParamFast (SCIP* scip, DEC_DETECTOR* detector, SCIP_RESULT* result):
    py_detector = get_py_detector(detector)
    result_dict = py_detector.setParamFast() or {}
    result[0] = result_dict.get("result", <SCIP_RESULT>result[0])
    return SCIP_OKAY
