cdef class Score:
    """Base class of the Score Plugin"""

    cdef public SCIPModel model
    cdef public str scorename

    def scorefree(self):
        '''calls destructor and frees memory of score'''
        pass

    def scorecalculate(self, partialdecid):
        '''calls calculate method of score'''
        return {}

cdef SCIP_RETCODE PyScoreFree(SCIP* scip, DEC_SCORE* score) with gil:
    cdef DEC_SCOREDATA* scoredata
    scoredata = GCGscoreGetData(score)
    py_score = <Score>scoredata
    py_score.scorefree()
    Py_DECREF(py_score)
    return SCIP_OKAY

cdef SCIP_RETCODE PyScoreCalculate(SCIP* scip, DEC_SCORE* score, int partialdecid, SCIP_Real* scorevalue) with gil:
    cdef DEC_SCOREDATA* scoredata
    scoredata = GCGscoreGetData(score)
    py_score = <Score>scoredata
    result_dict = py_score.scorecalculate(partialdecid)
    scorevalue[0] = result_dict.get("scorevalue", scorevalue[0])
    return SCIP_OKAY
