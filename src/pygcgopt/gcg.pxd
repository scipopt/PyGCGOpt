from pyscipopt.scip cimport SCIP, SCIP_RETCODE, SCIP_RESULT, SCIP_Bool, SCIP_Real, FILE, SCIP_CONS, SCIP_VAR, SCIP_PARAMSETTING

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp.map cimport map


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

    ctypedef struct DEC_DETECTORDATA:
        pass

    ctypedef struct DEC_DETECTOR:
        pass

    ctypedef struct PARTIALDEC_DETECTION_DATA:
        DETPROBDATA* detprobdata
        PARTIALDECOMP* workonpartialdec
        PARTIALDECOMP** newpartialdecs
        int nnewpartialdecs
        double detectiontime

    SCIP_RETCODE DECincludeDetector( SCIP* scip, const char* name, const char decchar, const char* description, int freqCallRound, int maxCallRound, int minCallRound, int freqCallRoundOriginal, int maxCallRoundOriginal, int minCallRoundOriginal, int priority, SCIP_Bool enabled, SCIP_Bool enabledFinishing, SCIP_Bool enabledPostprocessing, SCIP_Bool skip, SCIP_Bool usefulRecall, DEC_DETECTORDATA *detectordata, SCIP_RETCODE (*freeDetector) (SCIP* scip, DEC_DETECTOR* detector), SCIP_RETCODE (*initDetector) (SCIP* scip, DEC_DETECTOR* detector), SCIP_RETCODE (*exitDetector) (SCIP* scip, DEC_DETECTOR* detector), SCIP_RETCODE (*propagatePartialdecDetector) (SCIP* scip, DEC_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result), SCIP_RETCODE (*finishPartialdecDetector) (SCIP* scip, DEC_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result), SCIP_RETCODE (*postprocessPartialdecDetector) (SCIP* scip, DEC_DETECTOR* detector, PARTIALDEC_DETECTION_DATA* partialdecdetectiondata, SCIP_RESULT* result), SCIP_RETCODE (*setParamAggressiveDetector) (SCIP* scip, DEC_DETECTOR* detector, SCIP_RESULT* result), SCIP_RETCODE (*setParamDefaultDetector) (SCIP* scip, DEC_DETECTOR* detector, SCIP_RESULT* result), SCIP_RETCODE (*setParamFastDetector) (SCIP* scip, DEC_DETECTOR* detector, SCIP_RESULT* result))

    DEC_DETECTORDATA* DECdetectorGetData(DEC_DETECTOR* detector)

    DEC_DETECTOR** GCGconshdlrDecompGetDetectors(SCIP* scip)
    int GCGconshdlrDecompGetNDetectors(SCIP* scip)
    const char* DECdetectorGetName(DEC_DETECTOR* detector)

    SCIP* GCGgetMasterprob(SCIP* scip)

    ctypedef enum GP_OUTPUT_FORMAT:
        GP_OUTPUT_FORMAT_PDF
        GP_OUTPUT_FORMAT_PNG
        GP_OUTPUT_FORMAT_SVG


cdef extern from "gcg/pub_gcgsepa.h":
    SCIP_RETCODE GCGsetSeparators(SCIP* scip, SCIP_PARAMSETTING paramsetting)


cdef extern from "gcg/gcgplugins.h":
    SCIP_RETCODE SCIPincludeGcgPlugins(SCIP* scip)


cdef extern from "gcg/cons_decomp.h":
    SCIP_RETCODE GCGconshdlrDecompTranslateOrigPartialdecs(SCIP* scip)
    unsigned int GCGconshdlrDecompGetNPartialdecs(SCIP* scip)
    SCIP_RETCODE DECdetectStructure(SCIP* scip, SCIP_RESULT* result)
    SCIP_RETCODE SCIPconshdlrDecompRepairConsNames(SCIP* scip)
    unsigned int GCGconshdlrDecompGetNPartialdecs(SCIP* scip)
    SCIP_RETCODE GCGconshdlrDecompGetPartialdecsList(SCIP* scip, int** idlist, int* listlength)
    unsigned int GCGconshdlrDecompGetNDecomps(SCIP* scip)
    SCIP_RETCODE GCGconshdlrDecompGetFinishedPartialdecsList(SCIP* scip, int** idlist, int* listlength)
    SCIP_RETCODE DECwriteAllDecomps(SCIP* scip, char* directory, char* extension, SCIP_Bool original, SCIP_Bool presolved)


cdef extern from "gcg/cons_decomp.hpp":
    PARTIALDECOMP* GCGconshdlrDecompGetPartialdecFromID(SCIP* scip, int partialdecid)
    SCIP_RETCODE GCGconshdlrDecompAddPreexisitingPartialDec(SCIP* scip, PARTIALDECOMP* partialdec)


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
        vector[DEC_DETECTOR *] getDetectorchain() except +
        bool getFinishedByFinisher() except +
        unsigned long getHashValue() except +
        int getID() except +
        vector[int] getLinkingvars() except +
        vector[int] getMasterconss() except +
        vector[int] getMastervars() except +
        int getNCoeffsForBlock(int blockid) except +
        int getNCoeffsForMaster() except +
        # double getScore(SCORETYPE type) except +
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
        vector[int] getVarsForBlock(int block) except +
        int getVarProbindexForBlock(int varid, int block) except +
        bool isComplete() except +
        bool isConsMastercons(int cons) except +
        bool isConsOpencons(int cons) except +
        bool isAssignedToOrigProb() except +
        bool isSelected() except +
        SCIP_RETCODE isEqual(PARTIALDECOMP * otherpartialdec, unsigned int * isequal, bool sortpartialdecs) except +
        bool isPropagatedBy(DEC_DETECTOR * detector) except +
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
        void setDetectorchain(vector[DEC_DETECTOR *] givenDetectorChain) except +
        void setDetectorPropagated(DEC_DETECTOR * detector) except +
        void setDetectorFinished(DEC_DETECTOR * detector) except +
        void setDetectorFinishedOrig(DEC_DETECTOR * detectorID) except +
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
        void showVisualisation() except +
        void generateVisualisation(char * filename, char * outname, GP_OUTPUT_FORMAT outputformat) except +
        void writeVisualisationFile(char * filename, char * outname, GP_OUTPUT_FORMAT outputformat) except +
        map[pair[int, int], double] writeNonzerosWithRhsAndObj() except +
        unsigned int shouldCompletedByConsToMaster() except +
        bool sort() except +
        void setPctConssToBlockVector(vector[double] newvector) except +
        void setPctConssFromFreeVector(vector[double] newvector) except +
        void setPctConssToBorderVector(vector[double] newvector) except +
        void setPctVarsToBorderVector(vector[double] newvector) except +
        void setPctVarsToBlockVector(vector[double] newvector) except +
        void setPctVarsFromFreeVector(vector[double] newvector) except +
        void setDetectorClockTimes(vector[double] newvector) except +
        double getClassicScore() except +
        void setClassicScore(double score) except +
        double getBorderAreaScore() except +
        void setBorderAreaScore(double score) except +
        double getMaxWhiteScore() except +
        void setMaxWhiteScore(double score) except +
        double getMaxForWhiteScore() except +
        void setMaxForWhiteScore(double score) except +
        double getSetPartForWhiteScore() except +
        void setSetPartForWhiteScore(double score) except +
        double getMaxForWhiteAggScore() except +
        void setMaxForWhiteAggScore(double score) except +
        double getSetPartForWhiteAggScore() except +
        void setSetPartForWhiteAggScore(double score) except +
        double getBendersScore() except +
        void setBendersScore(double score) except +
        double getStrongDecompScore() except +
        void setStrongDecompScore(double score) except +
        void prepare() except +
        bool aggInfoCalculated() except +
        void calcAggregationInformation(bool ignoreDetectionLimits) except +
        vector[vector[int] ] getConssForBlocks() except +
        int getTranslatedpartialdecid() except +
        void setTranslatedpartialdecid(int decid) except +
        void buildDecChainString(char * buffer) except +

        bool fixConsToBlock(SCIP_CONS* cons, int block)
        bool fixConsToMaster(SCIP_CONS* cons)


cdef extern from "gcg/class_detprobdata.h" namespace "gcg":
    cdef cppclass DETPROBDATA:
        vector[pair[int, int] ] candidatesNBlocks
        vector[ConsPartition *] conspartitioncollection
        vector[VarPartition *] varpartitioncollection
        double classificationtime
        double nblockscandidatescalctime
        double postprocessingtime
        double translatingtime
        void addConsPartition(ConsPartition* partition) except +
        void addCandidatesNBlocksNVotes(int candidate, int nvotes) except +
        void addPartialdecToAncestor(PARTIALDECOMP * partialdec) except +
        bool addPartialdecToOpen(PARTIALDECOMP * partialdec) except +
        bool addPartialdecToFinished(PARTIALDECOMP * partialdec) except +
        void addPartialdecToFinishedUnchecked(PARTIALDECOMP * partialdec) except +
        void addVarPartition(VarPartition* partition) except +
        void clearAncestorPartialdecs() except +
        void clearCurrentPartialdecs() except +
        void clearFinishedPartialdecs() except +
        void createConssAdjacency() except +
        void freeTemporaryData() except +
        PARTIALDECOMP * getAncestorPartialdec(int partialdecindex) except +
        ConsPartition* getConsPartition(int partitionIndex) except +
        SCIP_CONS* getCons(int consIndex) except +
        vector[int] getConssForCons(int consIndex) except +
        vector[int] getConssForVar(int varIndex) except +
        vector[PARTIALDECOMP *] getOpenPartialdecs() except +
        PARTIALDECOMP * getFinishedPartialdec(int partialdecindex) except +
        vector[PARTIALDECOMP *] getFinishedPartialdecs() except +
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
        vector[SCIP_VAR *] getOrigVarsFixedZero() except +
        vector[SCIP_CONS *] getRelevantConss() except +
        vector[SCIP_VAR *] getRelevantVars() except +
        SCIP* getScip()
        void getSortedCandidatesNBlocks(vector[int] candidates) except +
        double getVal(int row, int col) except +
        vector[double] getValsForCons(int consIndex) except +
        VarPartition* getVarPartition(int partitionIndex) except +
        vector[VarPartition *] getVarPartitions() except +
        vector[int] getVarsForCons(int consIndex) except +
        bool isConsCardinalityCons(int consindexd) except +
        unsigned int isConssAdjInitialized() except +
        bool isConsSetpp(int consindexd) except +
        bool isConsSetppc(int consindexd) except +
        unsigned int isPartialdecDuplicateofFinished(PARTIALDECOMP * partialdec) except +
        unsigned int isAssignedToOrigProb() except +
        unsigned int partialdecIsNoDuplicateOfPartialdecs(PARTIALDECOMP * comppartialdec, vector[PARTIALDECOMP *] partialdecs, bool sort) except +
        void sortFinishedForScore() except +
        vector[PARTIALDECOMP *] translatePartialdecs(DETPROBDATA * otherdata, vector[PARTIALDECOMP *] otherpartialdecs) except +


cdef extern from "gcg/class_conspartition.h" namespace "gcg":
    cdef cppclass ConsPartition:
        ConsPartition(ConsPartition * toCopy)
        void assignConsToClass(int consindex, int classindex) except +
        vector[vector[int] ] getAllSubsets(bool both, bool only_master, bool only_pricing) except +
        char * getClassNameOfCons(int consindex)
        int getClassOfCons(int consindex) except +
        int getNConss() except +
        vector[int] getNConssOfClasses() except +
        bool isConsClassified(int consindex) except +
        ConsPartition * reduceClasses(int maxNumberOfClasses) except +
        const char* getName() except +


cdef extern from "gcg/class_varpartition.h" namespace "gcg":
    cdef cppclass VarPartition:
        VarPartition(VarPartition * toCopy)
        void assignVarToClass(int varindex, int classindex) except +
        vector[vector[int] ] getAllSubsets(bool all, bool linking, bool master, bool block) except +
        char * getClassNameOfVar(int varindex)
        int getClassOfVar(int varindex) except +
        int getNVars() except +
        vector[int] getNVarsOfClasses() except +
        bool isVarClassified(int varindex) except +
        VarPartition * reduceClasses(int maxNumberOfClasses) except +


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
