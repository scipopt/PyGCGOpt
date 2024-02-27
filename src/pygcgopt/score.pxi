cdef class Score:
    """Base class of the Score Plugin"""

    cdef public Model model
    cdef public str scorename

    def scorefree(self):
        '''calls destructor and frees memory of score'''
        pass

    def scorecalculate(self, partialdec):
        '''calls calculate method of score'''
        return {}

cdef SCIP_RETCODE PyScoreFree(SCIP* scip, GCG_SCORE* score) noexcept with gil:
    cdef GCG_SCOREDATA* scoredata
    scoredata = GCGscoreGetData(score)
    py_score = <Score>scoredata
    py_score.scorefree()
    Py_DECREF(py_score)
    return SCIP_OKAY

cdef SCIP_RETCODE PyScoreCalculate(SCIP* scip, GCG_SCORE* score, int partialdecid, SCIP_Real* scorevalue) noexcept with gil:
    cdef GCG_SCOREDATA* scoredata
    scoredata = GCGscoreGetData(score)
    py_score = <Score>scoredata
    partialdec = py_score.model.getPartDecompFromID(partialdecid)
    result_dict = py_score.scorecalculate(partialdec)
    scorevalue[0] = result_dict.get("scorevalue", scorevalue[0])
    return SCIP_OKAY
