from pyscipopt.scip cimport SCIP, SCIP_RETCODE, SCIP_RESULT, SCIP_Bool, SCIP_Real, FILE, SCIP_CONS, SCIP_VAR, SCIP_PARAMSETTING, SCIP_SOL

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair


cdef extern from "limits.h":
    cdef int INT_MAX


cdef extern from "gcg/gcg.h":
    void GCGprintVersion(SCIP* scip, FILE* file)
    SCIP_RETCODE GCGprintStatistics(SCIP* scip, FILE* file)
    SCIP_RETCODE GCGtransformProb(SCIP* scip)
    SCIP_RETCODE GCGpresolve(SCIP* scip)
    SCIP_RETCODE GCGdetect(SCIP* scip)
    SCIP_RETCODE GCGsolve(SCIP* scip)
    SCIP_Bool GCGdetectionTookPlace(SCIP* scip, SCIP_Bool original)

    SCIP_Real GCGgetDualbound(SCIP* scip)

    ctypedef struct GCG_DETECTORDATA:
        pass

    ctypedef struct GCG_DETECTOR:
        pass

    ctypedef struct GCG_SCOREDATA:
        pass

    ctypedef struct GCG_SCORE:
        pass

    ctypedef struct PARTIALDEC_DETECTION_DATA:
        DETPROBDATA* detprobdata
        PARTIALDECOMP* workonpartialdec
        PARTIALDECOMP** newpartialdecs
        int nnewpartialdecs
        double detectiontime

    SCIP_RETCODE GCGincludeDetector(SCIP* scip, const char* name, const char decchar, const char* description, int freqCallRound, int maxCallRound, int minCallRound, int freqCallRoundOriginal, int maxCallRoundOriginal, int minCallRoundOriginal, int priority, SCIP_Bool enabled, SCIP_Bool enabledFinishing, SCIP_Bool enabledPostprocessing, SCIP_Bool skip, SCIP_Bool usefulRecall, GCG_DETECTORDATA *detectordata, SCIP_RETCODE (*freeDetector) (SCIP* scip, GCG_DETECTOR* detector), SCIP_RETCODE (*initDetector) (SCIP* scip, GCG_DETECTOR* detector), SCIP_RETCODE (*exitDetector) (SCIP* scip, GCG_DETECTOR* detector), SCIP_RETCODE (*propagatePartialdecDetector) (SCIP* scip, GCG_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result), SCIP_RETCODE (*finishPartialdecDetector) (SCIP* scip, GCG_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result), SCIP_RETCODE (*postprocessPartialdecDetector) (SCIP* scip, GCG_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result), SCIP_RETCODE (*setParamAggressiveDetector) (SCIP* scip, GCG_DETECTOR* detector, SCIP_RESULT* result), SCIP_RETCODE (*setParamDefaultDetector) (SCIP* scip, GCG_DETECTOR* detector, SCIP_RESULT* result), SCIP_RETCODE (*setParamFastDetector) (SCIP* scip, GCG_DETECTOR* detector, SCIP_RESULT* result))

    GCG_DETECTORDATA* GCGdetectorGetData(GCG_DETECTOR* detector)

    GCG_DETECTOR** GCGconshdlrDecompGetDetectors(SCIP* scip)
    int GCGconshdlrDecompGetNDetectors(SCIP* scip)
    const char* GCGdetectorGetName(GCG_DETECTOR* detector)

    SCIP* GCGgetMasterprob(SCIP* scip)

    #ConsClassifier
    ctypedef struct GCG_CLASSIFIERDATA:
        pass

    ctypedef struct GCG_CONSCLASSIFIER:
        pass

    SCIP_RETCODE GCGincludeConsClassifier(SCIP* scip, const char* name, const char* description, int priority, SCIP_Bool enabled, GCG_CLASSIFIERDATA *classifierdata, SCIP_RETCODE (*freeClassifier) (SCIP* scip, GCG_CONSCLASSIFIER* classifier), SCIP_RETCODE (*classify) (SCIP* scip, GCG_CONSCLASSIFIER* classifierpointer, SCIP_Bool transformed))

    GCG_CLASSIFIERDATA* GCGconsClassifierGetData(GCG_CONSCLASSIFIER* classifier)

    GCG_CONSCLASSIFIER** GCGconshdlrDecompGetConsClassifiers(SCIP* scip)
    int GCGconshdlrDecompGetNConsClassifiers(SCIP* scip)
    const char* GCGconsClassifierGetName(GCG_CONSCLASSIFIER* classifier)

    #VarClassifier
    ctypedef struct GCG_VARCLASSIFIER:
        pass

    SCIP_RETCODE GCGincludeVarClassifier(SCIP* scip, const char* name, const char* description, int priority, SCIP_Bool enabled, GCG_CLASSIFIERDATA* classifierdata, SCIP_RETCODE (*freeClassifier) (SCIP* scip, GCG_VARCLASSIFIER* classifier), SCIP_RETCODE (*classify) (SCIP* scip, GCG_VARCLASSIFIER* classifierpointer, SCIP_Bool transformed))

    GCG_CLASSIFIERDATA* GCGvarClassifierGetData(GCG_VARCLASSIFIER* classifier)

    int GCGconshdlrDecompGetNVarClassifiers(SCIP* scip)
    GCG_VARCLASSIFIER** GCGconshdlrDecompGetVarClassifiers(SCIP* scip)
    const char* GCGvarClassifierGetName(GCG_VARCLASSIFIER* classifier)

    ctypedef enum GP_OUTPUT_FORMAT:
        GP_OUTPUT_FORMAT_PDF
        GP_OUTPUT_FORMAT_PNG
        GP_OUTPUT_FORMAT_SVG

    SCIP_RETCODE GCGincludeScore(SCIP* scip, const char* name, const char* shortname,const char* description, GCG_SCOREDATA* scoredata, SCIP_RETCODE (*scorefree) (SCIP* scip, GCG_SCORE* score), SCIP_RETCODE (*scorecalc) (SCIP* scip, GCG_SCORE* score, int partialdecid, SCIP_Real* scorevalue))
    GCG_SCOREDATA* GCGscoreGetData(GCG_SCORE* score)
    int GCGgetNScores(SCIP* scip)
    GCG_SCORE** GCGgetScores(SCIP* scip)
    const char* GCGscoreGetName(GCG_SCORE* score)
    GCG_SCORE* GCGfindScore(SCIP* scip, const char* name)


cdef extern from "gcg/pub_gcgsepa.h":
    SCIP_RETCODE GCGsetSeparators(SCIP* scip, SCIP_PARAMSETTING paramsetting)


cdef extern from "gcg/gcgplugins.h":
    SCIP_RETCODE SCIPincludeGcgPlugins(SCIP* scip)

cdef extern from "gcg/pricer_gcg.h":
    ctypedef struct GCG_SOLVER:
        pass

    ctypedef struct GCG_SOLVERDATA:
        pass

    ctypedef enum GCG_PRICINGSTATUS:
        GCG_PRICINGSTATUS_UNKNOWN        = 0
        GCG_PRICINGSTATUS_NOTAPPLICABLE  = 1
        GCG_PRICINGSTATUS_SOLVERLIMIT    = 2
        GCG_PRICINGSTATUS_OPTIMAL        = 3
        GCG_PRICINGSTATUS_INFEASIBLE     = 4
        GCG_PRICINGSTATUS_UNBOUNDED      = 5

    SCIP_RETCODE GCGpricerIncludeSolver(
        SCIP* scip,
        const char* name,
        const char* desc,
        int priority,
        SCIP_Bool heurenabled,
        SCIP_Bool exactenabled,
        SCIP_RETCODE (*solverupdate) (SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Bool varobjschanged, SCIP_Bool varbndschanged, SCIP_Bool consschanged),
        SCIP_RETCODE (*solversolve) (SCIP* scip, SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Real dualsolconv, SCIP_Real* lowerbound, GCG_PRICINGSTATUS* status),
        SCIP_RETCODE (*solveheur) (SCIP* scip, SCIP* pricingprob, GCG_SOLVER* solver, int probnr, SCIP_Real dualsolconv, SCIP_Real* lowerbound, GCG_PRICINGSTATUS* status),
        SCIP_RETCODE (*solverfree) (SCIP* scip, GCG_SOLVER* solver),
        SCIP_RETCODE (*solverinit) (SCIP* scip, GCG_SOLVER* solver),
        SCIP_RETCODE (*solverexit) (SCIP* scip, GCG_SOLVER* solver),
        SCIP_RETCODE (*solverinitsol) (SCIP* scip, GCG_SOLVER* solver),
        SCIP_RETCODE (*solverexitsol) (SCIP* scip, GCG_SOLVER* solver),
        GCG_SOLVERDATA*       solverdata
    )
    GCG_SOLVER** GCGpricerGetSolvers(SCIP* scip)
    int GCGpricerGetNSolvers(SCIP* scip)
    SCIP_RETCODE GCGpricerAddCol(SCIP* scip, GCG_COL* col)
    const char* GCGsolverGetName(GCG_SOLVER* solver)


cdef extern from "gcg/pub_solver.h":
    GCG_SOLVERDATA* GCGsolverGetData(GCG_SOLVER* solver)


cdef extern from "gcg/pub_gcgcol.h":
    SCIP_RETCODE GCGcreateGcgCol(SCIP* scip, GCG_COL** gcgcol, int prob, SCIP_VAR** vars, SCIP_Real* vals, int nvars, SCIP_Bool isray, SCIP_Real redcost)
    SCIP_RETCODE GCGcreateGcgColFromSol(SCIP* scip, GCG_COL** gcgcol, int prob, SCIP_SOL* sol, SCIP_Bool isray, SCIP_Real redcost)

    ctypedef struct GCG_COL:
        pass


cdef extern from "gcg/cons_decomp.h":
    SCIP_RETCODE GCGconshdlrDecompTranslateOrigPartialdecs(SCIP* scip)
    unsigned int GCGconshdlrDecompGetNPartialdecs(SCIP* scip)
    SCIP_RETCODE DECdetectStructure(SCIP* scip, SCIP_RESULT* result)
    SCIP_RETCODE SCIPconshdlrDecompRepairConsNames(SCIP* scip)
    unsigned int GCGconshdlrDecompGetNPartialdecs(SCIP* scip)
    SCIP_RETCODE GCGconshdlrDecompGetPartialdecsList(SCIP* scip, int** idlist, int* listlength)
    unsigned int GCGconshdlrDecompGetNDecomps(SCIP* scip)
    SCIP_RETCODE GCGconshdlrDecompGetFinishedPartialdecsList(SCIP* scip, int** idlist, int* listlength)
    SCIP_RETCODE GCGwriteAllDecomps(SCIP* scip, char* directory, char* extension, SCIP_Bool original, SCIP_Bool presolved)


cdef extern from "gcg/cons_decomp.hpp":
    PARTIALDECOMP* GCGconshdlrDecompGetPartialdecFromID(SCIP* scip, int partialdecid)
    SCIP_RETCODE GCGconshdlrDecompAddPreexisitingPartialDec(SCIP* scip, PARTIALDECOMP* partialdec)
    DETPROBDATA* GCGconshdlrDecompGetDetprobdataPresolved(SCIP* scip)
    DETPROBDATA* GCGconshdlrDecompGetDetprobdataOrig(SCIP* scip)


cdef extern from "gcg/class_partialdecomp.h" namespace "gcg":
    ctypedef enum USERGIVEN:
        NOT
        PARTIAL
        COMPLETE
        COMPLETED_CONSTOMASTER

    cdef cppclass PARTIALDECOMP:
        PARTIALDECOMP(SCIP * scip, bool originalProblem)
        PARTIALDECOMP(PARTIALDECOMP * partialdecToCopy)
        int addBlock() except +
        void addClockTime(double clocktime) except +
        void addDecChangesFromAncestor(PARTIALDECOMP * ancestor) except +
        void addDetectorChainInfo(char * decinfo) except +
        void addNNewBlocks(int nnewblocks) except +
        void addPctConssFromFree(double pct) except +
        void addPctConssToBlock(double pct) except +
        void addPctConssToBorder(double pct) except +
        void addPctVarsFromFree(double pct) except +
        void addPctVarsToBlock(double pct) except +
        void addPctVarsToBorder(double pct) except +
        bool alreadyAssignedConssToBlocks() except +
        # SCIP_RETCODE assignBorderFromConstoblock(SCIP_HASHMAP * constoblock, int givenNBlocks) except +
        bool assignCurrentStairlinking() except +
        void assignOpenConssToMaster() except +
        # SCIP_RETCODE assignPartialdecFromConstoblock(SCIP_HASHMAP * constoblock, int additionalNBlocks) except +
        SCIP_RETCODE assignPartialdecFromConstoblockVector(vector[int] constoblock, int additionalNBlocks) except +
        void assignSmallestComponentsButOneConssAdjacency() except +
        void calcStairlinkingVars() except +
        bool checkAllConssAssigned() except +
        bool checkConsistency() except +
        void complete() except +
        void completeByConnected() except +
        void completeByConnectedConssAdjacency() except +
        void completeGreedily() except +
        void removeMastercons(int consid) except +
        void considerImplicits() except +
        void copyPartitionStatistics(PARTIALDECOMP * otherpartialdec) except +
        void deleteEmptyBlocks(bool variables) except +
        void deleteOpencons(int opencons) except +
        void deleteOpenvar(int openvar) except +
        void displayInfo(int detailLevel) except +
        # SCIP_RETCODE filloutBorderFromConstoblock(SCIP_HASHMAP * constoblock, int givenNBlocks) except +
        # SCIP_RETCODE filloutPartialdecFromConstoblock(SCIP_HASHMAP * constoblock, int givenNBlocks) except +
        void findVarsLinkingToMaster() except +
        void findVarsLinkingToStairlinking() except +
        int getAncestorID(int ancestorindex) except +
        vector[int] getAncestorList() except +
        void setAncestorList(vector[int] newlist) except +
        void removeAncestorID(int ancestorid) except +
        void addAncestorID(int ancestor) except +
        vector[int] getBlocksForRep(int repid) except +
        double getDetectorClockTime(int detectorchainindex) except +
        vector[double] getDetectorClockTimes() except +
        vector[int] getConssForBlock(int block) except +
        vector[GCG_DETECTOR *] getDetectorchain() except +
        bool getFinishedByFinisher() except +
        unsigned long getHashValue() except +
        int getID() except +
        vector[int] getLinkingvars() except +
        vector[int] getMasterconss() except +
        vector[int] getMastervars() except +
        int getNCoeffsForBlock(int blockid) except +
        int getNCoeffsForMaster() except +
        unsigned int hasSetppccardMaster() except +
        unsigned int hasSetppcMaster() except +
        unsigned int hasSetppMaster() except +
        USERGIVEN getUsergiven() except +
        int getNAncestors() except +
        int getNBlocks() except +
        int getNConss() except +
        int getNConssForBlock(int block) except +
        vector[string] getDetectorchainInfo() except +
        int getNDetectors() except +
        int getNLinkingvars() except +
        int getNMasterconss() except +
        int getNMastervars() except +
        int getNNewBlocks(int detectorchainindex) except +
        vector[int] getNNewBlocksVector() except +
        int getNTotalStairlinkingvars() except +
        int getNOpenconss() except +
        int getNOpenvars() except +
        int getNReps() except +
        int getNStairlinkingvars(int block) except +
        int getNVars() except +
        int getNVarsForBlock(int block) except +
        int getNVarsForBlocks() except +
        int * getOpenconss() except +
        vector[int] getOpenconssVec() except +
        int * getOpenvars() except +
        vector[int] getOpenvarsVec() except +
        double getPctVarsToBorder(int detectorchainindex) except +
        vector[double] getPctVarsToBorderVector() except +
        double getPctVarsToBlock(int detectorchainindex) except +
        vector[double] getPctVarsToBlockVector() except +
        double getPctVarsFromFree(int detectorchainindex) except +
        vector[double] getPctVarsFromFreeVector() except +
        double getPctConssToBorder(int detectorchainindex) except +
        vector[double] getPctConssToBorderVector() except +
        double getPctConssToBlock(int detectorchainindex) except +
        vector[double] getPctConssToBlockVector() except +
        double getPctConssFromFree(int detectorchainindex) except +
        vector[double] getPctConssFromFreeVector() except +
        int getRepForBlock(int blockid) except +
        vector[int] getRepVarmap(int repid, int blockrepid) except +
        DETPROBDATA * getDetprobdata() except +
        int * getStairlinkingvars(int block) except +
        #vector[int] getStairlinkingvarsVec(int block) except +
        vector[int] getVarsForBlock(int block) except +
        int getVarProbindexForBlock(int varid, int block) except +
        bool isComplete() except +
        bool isConsMastercons(int cons) except +
        bool isConsOpencons(int cons) except +
        bool isAssignedToOrigProb() except +
        bool isSelected() except +
        SCIP_RETCODE isEqual(PARTIALDECOMP * otherpartialdec, unsigned int * isequal, bool sortpartialdecs) except +
        bool isPropagatedBy(GCG_DETECTOR * detector) except +
        bool isTrivial() except +
        bool isVarBlockvarOfBlock(int var, int block) except +
        bool isVarLinkingvar(int var) except +
        bool isVarMastervar(int var) except +
        bool isVarOpenvar(int var) except +
        bool isVarStairlinkingvar(int var) except +
        bool isVarStairlinkingvarOfBlock(int var, int block) except +
        void printPartitionInformation(SCIP * givenscip, FILE * file) except +
        void refineToBlocks() except +
        void refineToMaster() except +
        void setConsPartitionStatistics(int detectorchainindex, ConsPartition * partition, vector[int] consclassesmaster) except +
        void setConsToBlock(int consToBlock, int block) except +
        # void fixConsToBlock(int cons, int block) except +
        void setConsToMaster(int consToMaster) except +
        # void fixConsToMaster(int cons) except +
        void setDetectorchain(vector[GCG_DETECTOR *] givenDetectorChain) except +
        void setDetectorPropagated(GCG_DETECTOR * detector) except +
        void setDetectorFinished(GCG_DETECTOR * detector) except +
        void setDetectorFinishedOrig(GCG_DETECTOR * detectorID) except +
        void setFinishedByFinisher(bool finished) except +
        void setFinishedByFinisherOrig(bool finished) except +
        void setNBlocks(int nblocks) except +
        void setSelected(bool selected) except +
        void setStemsFromOrig(bool fromorig) except +
        void setUsergiven(USERGIVEN usergiven) except +
        void setVarPartitionStatistics(int detectorchainindex, VarPartition * partition, vector[int] varclasseslinking, vector[int] varclassesmaster) except +
        void setVarToBlock(int varToBlock, int block) except +
        void fixVarToBlock(int var, int block) except +
        void setVarToLinking(int varToLinking) except +
        void fixVarToLinking(int var) except +
        void setVarToMaster(int varToMaster) except +
        void fixVarToMaster(int var) except +
        void setVarToStairlinking(int varToStairLinking, int block1, int block2) except +
        void fixVarToStairlinking(int var, int firstblock) except +
        bool fixConsToBlockByName(char * consname, int blockid) except +
        bool fixVarToBlockByName(char * varname, int blockid) except +
        bool fixConsToMasterByName(char * consname) except +
        bool fixVarToMasterByName(char * varname) except +
        bool fixVarToLinkingByName(char * varname) except +
        void showVisualization() except +
        void generateVisualization(char * filename, char * outname, GP_OUTPUT_FORMAT outputformat) except +
        void writeVisualizationFile(char * filename, char * outname, GP_OUTPUT_FORMAT outputformat) except +
        unsigned int shouldCompletedByConsToMaster() except +
        bool sort() except +
        void setPctConssToBlockVector(vector[double] newvector) except +
        void setPctConssFromFreeVector(vector[double] newvector) except +
        void setPctConssToBorderVector(vector[double] newvector) except +
        void setPctVarsToBorderVector(vector[double] newvector) except +
        void setPctVarsToBlockVector(vector[double] newvector) except +
        void setPctVarsFromFreeVector(vector[double] newvector) except +
        void setDetectorClockTimes(vector[double] newvector) except +
        double getMaxWhiteScore() except +
        void prepare() except +
        bool aggInfoCalculated() except +
        void calcAggregationInformation(bool ignoreDetectionLimits) except +
        vector[vector[int]] getConssForBlocks() except +
        int getTranslatedpartialdecid() except +
        void setTranslatedpartialdecid(int decid) except +
        void buildDecChainString(char * buffer) except +
        bool fixConsToBlock(SCIP_CONS* cons, int block)
        bool fixConsToMaster(SCIP_CONS* cons)
        SCIP_Real getScore(GCG_SCORE* score) except +


cdef extern from "gcg/class_detprobdata.h" namespace "gcg":
    cdef cppclass DETPROBDATA:
        vector[pair[int, int]] candidatesNBlocks
        vector[ConsPartition *] conspartitioncollection
        vector[VarPartition *] varpartitioncollection
        double classificationtime
        double nblockscandidatescalctime
        double postprocessingtime
        double translatingtime
        DETPROBDATA(SCIP* scip, SCIP_Bool _originalProblem)
        void addConsPartition(ConsPartition* partition) except +
        void addCandidatesNBlocksNVotes(int candidate, int nvotes) except +
        void addPartialdecToAncestor(PARTIALDECOMP* partialdec) except +
        bool addPartialdecToOpen(PARTIALDECOMP* partialdec) except +
        bool addPartialdecToFinished(PARTIALDECOMP* partialdec) except +
        void addPartialdecToFinishedUnchecked(PARTIALDECOMP* partialdec) except +
        void addVarPartition(VarPartition* partition) except +
        void clearAncestorPartialdecs() except +
        void clearCurrentPartialdecs() except +
        void clearFinishedPartialdecs() except +
        void createConssAdjacency() except +
        void freeTemporaryData() except +
        PARTIALDECOMP* getAncestorPartialdec(int partialdecindex) except +
        ConsPartition* getConsPartition(int partitionIndex) except +
        SCIP_CONS* getCons(int consIndex) except +
        vector[int] getConssForCons(int consIndex) except +
        vector[int] getConssForVar(int varIndex) except +
        vector[PARTIALDECOMP*] getOpenPartialdecs() except +
        PARTIALDECOMP* getFinishedPartialdec(int partialdecindex) except +
        vector[PARTIALDECOMP*] getFinishedPartialdecs() except +
        int getIndexForCons(SCIP_CONS* cons) except +
        int getIndexForVar(SCIP_VAR* var) except +
        int getNAncestorPartialdecs() except +
        int getNConsPartitions() except +
        int getNConss() except +
        int getNConssForCons(int consIndex) except +
        int getNConssForVar(int varIndex) except +
        int getNOpenPartialdecs() except +
        int getNFinishedPartialdecs() except +
        int getNPartialdecs() except +
        int getNNonzeros() except +
        int getNVarPartitions() except +
        int getNVars() except +
        int getNVarsForCons(int consIndex) except +
        vector[SCIP_VAR*] getOrigVarsFixedZero() except +
        vector[SCIP_CONS*] getRelevantConss() except +
        vector[SCIP_VAR*] getRelevantVars() except +
        SCIP* getScip()
        void getSortedCandidatesNBlocks(vector[int] candidates) except +
        SCIP_Real getVal(int row, int col) except +
        vector[SCIP_Real] getValsForCons(int consIndex) except +
        VarPartition* getVarPartition(int partitionIndex) except +
        vector[VarPartition*] getVarPartitions() except +
        SCIP_VAR* getVar(int varIndex) except +
        vector[int] getVarsForCons(int consIndex) except +
        bool isConsCardinalityCons(int consindexd) except +
        SCIP_Bool isConssAdjInitialized() except +
        bool isConsSetpp(int consindexd) except +
        bool isConsSetppc(int consindexd) except +
        SCIP_Bool isPartialdecDuplicateofFinished(PARTIALDECOMP * partialdec) except +
        SCIP_Bool isAssignedToOrigProb() except +
        SCIP_Bool partialdecIsNoDuplicateOfPartialdecs(PARTIALDECOMP* comppartialdec, vector[PARTIALDECOMP*] partialdecs, bool sort) except +
        void sortFinishedForScore() except +
        vector[PARTIALDECOMP*] translatePartialdecs(DETPROBDATA* otherdata, vector[PARTIALDECOMP*] otherpartialdecs) except +


cdef extern from "gcg/class_conspartition.h" namespace "gcg":
    ctypedef enum CONS_DECOMPINFO:
        BOTH         = 0
        ONLY_MASTER  = 1
        ONLY_PRICING = 2

    cdef cppclass ConsPartition:
        # methods of superclass IndexPartition
        const char* getClassDescription(int classindex) except +
        const char* getClassName(int classindex) except +
        const char* getName() except +
        int getNClasses() except +
        #vector[int] reduceClasses(int maxNumberOfClasses) except +
        int removeEmptyClasses() except +
        void setClassDescription(int classindex, const char* desc) except +
        void setClassName(int classindex, const char* name) except +

        # constructors and methods of class ConsPartition
        ConsPartition(SCIP* scip, const char* name, int nClasses, int nConss) except +
        ConsPartition(ConsPartition * toCopy) except +
        int addClass(const char* name, const char* desc, CONS_DECOMPINFO decompInfo) except +
        void assignConsToClass(int consindex, int classindex) except +
        vector[vector[int]] getAllSubsets(bool both, bool only_master, bool only_pricing) except +
        CONS_DECOMPINFO getClassDecompInfo(int classindex) except +
        const char* getClassNameOfCons(int consindex) except +
        int getClassOfCons(int consindex) except +
        #const int* getConssToClasses() except +
        int getNConss() except +
        vector[int] getNConssOfClasses() except +
        bool isConsClassified(int consindex) except +
        ConsPartition* reduceClasses(int maxNumberOfClasses) except +
        void setClassDecompInfo(int classindex, CONS_DECOMPINFO decompInfo) except +

cdef extern from "gcg/class_varpartition.h" namespace "gcg":
    ctypedef enum VAR_DECOMPINFO:
        ALL     = 0
        LINKING = 1
        MASTER  = 2
        BLOCK   = 3

    cdef cppclass VarPartition:
        # methods of superclass IndexPartition
        const char* getClassDescription(int classindex) except +
        const char* getClassName(int classindex) except +
        const char* getName() except +
        int getNClasses() except +
        #vector[int] reduceClasses(int maxNumberOfClasses) except +
        int removeEmptyClasses() except +
        void setClassDescription(int classindex, const char* desc) except +
        void setClassName(int classindex, const char* name) except +

        # constructors and methods of class VarPartition
        VarPartition(SCIP* scip, const char* name, int nClasses, int nVars) except +
        VarPartition(VarPartition * toCopy) except +
        int addClass(const char* name, const char* desc, VAR_DECOMPINFO decompInfo) except +
        void assignVarToClass(int varindex, int classindex) except +
        vector[vector[int]] getAllSubsets(bool all, bool linking, bool master, bool block) except +
        VAR_DECOMPINFO getClassDecompInfo(int classindex) except +
        const char* getClassNameOfVar(int varindex) except +
        int getClassOfVar(int varindex) except +
        #const int* getVarsToClasses() except +
        int getNVars() except +
        vector[int] getNVarsOfClasses() except +
        bool isVarClassified(int varindex) except +
        VarPartition* reduceClasses(int maxNumberOfClasses) except +
        void setClassDecompInfo(int classindex, VAR_DECOMPINFO decompInfo) except +

cdef extern from "scip/scip.h":
    ctypedef struct SCIP_CLOCK:
        pass

    # Timing Functions
    SCIP_RETCODE SCIPcreateClock(SCIP* scip, SCIP_CLOCK** clck)
    SCIP_RETCODE SCIPfreeClock(SCIP* scip, SCIP_CLOCK** clck)
    SCIP_RETCODE SCIPresetClock(SCIP* scip, SCIP_CLOCK* clck)
    SCIP_RETCODE SCIPstartClock(SCIP* scip, SCIP_CLOCK* clck)
    SCIP_RETCODE SCIPstopClock(SCIP* scip, SCIP_CLOCK* clck)
    SCIP_Real SCIPgetClockTime(SCIP* scip, SCIP_CLOCK* clck)

cdef class DetProbData:
    cdef DETPROBDATA* thisptr

    @staticmethod
    cdef create(DETPROBDATA* thisptr)

cdef class ConsPart:
    cdef ConsPartition* consPartition
    cdef DetProbData detProbData

    @staticmethod
    cdef create(ConsPartition* thisptr, DetProbData detProbData)

cdef class VarPart:
    cdef VarPartition* varPartition
    cdef DetProbData detProbData

    @staticmethod
    cdef create(VarPartition* thisptr, DetProbData detProbData)

cdef class GCGColumn:
    cdef GCG_COL* gcg_col

    @staticmethod
    cdef create(GCG_COL* gcgcol)
