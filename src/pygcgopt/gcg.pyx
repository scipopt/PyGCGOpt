# distutils: language = c++

from pyscipopt.scip import PY_SCIP_CALL
from pyscipopt.scip cimport Model, Variable, Constraint, SCIP_RESULT, SCIP_DIDNOTRUN, SCIPgetStage, SCIP_STAGE, SCIP_STAGE_PRESOLVED, SCIP_OKAY

from cpython cimport Py_INCREF, Py_DECREF

from libc.stdlib cimport malloc, free

from libcpp cimport bool

from typing import List

from pygcgopt.util import str_conversion

from pathlib import Path
import tempfile
import weakref
from copy import copy

from collections.abc import Iterable

include "detector.pxi"


cdef SCIP_CLOCK* start_new_clock(SCIP* scip):
    cdef SCIP_CLOCK* clock
    PY_SCIP_CALL( SCIPcreateClock(scip, &clock) )
    PY_SCIP_CALL( SCIPstartClock(scip, clock) )
    return clock


cdef double stop_and_free_clock(SCIP* scip, SCIP_CLOCK* clock):
    PY_SCIP_CALL( SCIPstopClock(scip, clock) )
    cdef double detection_time = SCIPgetClockTime(scip, clock)
    PY_SCIP_CALL(SCIPfreeClock(scip, &clock) )
    return detection_time


cdef class GCGModel(Model):
    """!Main class for interaction with the GCG solver."""

    def includeDefaultPlugins(self):
        """!@brief Includes all default plug-ins of GCG into SCIP
        
        Called automatically during initialization of the model.
        """
        PY_SCIP_CALL(SCIPincludeGcgPlugins(self._scip))

    def presolve(self):
        """!@brief Presolve the problem."""
        PY_SCIP_CALL( GCGpresolve(self._scip) )

    def detect(self):
        """!@brief Detect the problem.
        
        Can be executed before or after presolving. If executed before presolving, the structure is detected on the original problem and presolving is skiped when solving the problem later.

        @see presolve()
        @see optimize()
        """
        PY_SCIP_CALL( GCGdetect(self._scip) )

    def printStatistics(self):
        """!@brief Print solving statistics of GCG to stdout."""
        PY_SCIP_CALL( GCGprintStatistics(self._scip, NULL) )

    def printVersion(self):
        """Print version, copyright information and compile mode of GCG and SCIP"""
        GCGprintVersion(self._scip, NULL)

        super().printVersion()

    def optimize(self):
        """!@brief Optimize the problem.
        
        This will transform, presolve and detect the problem if neccessary.
        Otherwise, GCG will solve the problem directly."""
        PY_SCIP_CALL( GCGsolve(self._scip) )

    def getDualbound(self):
        """!@brief Retrieve the best dual bound.
        
        This retrieves the same dual bound that GCG reports in the console log. The dual bound is based on the
        objective value of the optimized linear programming relaxation at the current node.

        @note The dual bound at the root node is *not* always equal to the solution of the restricted master problem LP
        relaxation. This can be due to master cuts or abortion of the pricing loop *before* the restricted master
        problem is optimal.

        @return The best dual bound of the current node.
        """
        return GCGgetDualbound(self._scip)

    def listDecompositions(self) -> List[PartialDecomposition]:
        """!@brief Lists all finnished decompositions found during the detection loop or provided by the user."""
        cdef int npartialdecs = GCGconshdlrDecompGetNPartialdecs(self._scip)
        cdef int* decids = <int*>malloc(npartialdecs * sizeof(int))

        GCGconshdlrDecompGetFinishedPartialdecsList(self._scip, &decids, &npartialdecs)

        decomps = [PartialDecomposition.create(GCGconshdlrDecompGetPartialdecFromID(self._scip, decids[i])) for i in range(npartialdecs)]

        free(decids)

        return decomps

    def addPreexistingPartialDecomposition(self, PartialDecomposition partialdec):
        """!@brief Adds a user specified decomposition to GCG.

        The passed PartialDecomposition can be partial or finnished. A partial decomposition will be completed by GCG using
        its detector loop. If a finnished decomposition is passed, GCG will skip the detection loop and use the
        provided decomposition right away.
        """
        partialdec.prepare()
        partialdec.setUsergiven()
        GCGconshdlrDecompAddPreexisitingPartialDec(self._scip, partialdec.thisptr)

    def createPartialDecomposition(self):
        """!@brief Creates a new empty PartialDecomposition.

        The created PartialDecomposition object can be used to fix constraints and variables. Afterwards, it can be
        passed to the model through addPreexistingPartialDecomposition().

        @see PartialDecomposition#fixConsToMaster()
        @see PartialDecomposition#fixConssToMaster()
        @see PartialDecomposition#fixConsToBlock()
        @see PartialDecomposition#fixConssToBlock()
        @see PartialDecomposition#fixConsToBlockId()
        @see PartialDecomposition#fixConssToBlockId()
        """
        cdef bool is_presolved = self.getStage() >= SCIP_STAGE_PRESOLVED
        cdef PARTIALDECOMP *decomp = new PARTIALDECOMP(self._scip, not is_presolved)
        return PartialDecomposition.create(decomp)

    def includeDetector(self, Detector detector, detectorname, decchar, desc, freqcallround=1, maxcallround=INT_MAX, mincallround=0, freqcallroundoriginal=1, maxcallroundoriginal=INT_MAX, mincallroundoriginal=0, priority=0, enabled=True, enabledfinishing=False, enabledpostprocessing=False, skip=False, usefulrecall=False):
        """!@brief includes a detector
        @param detector An object of a subclass of detector#Detector.
        @param detectorname name of the detector

        For an explanation for all arguments, see @ref DECincludeDetector().
        """
        if len(decchar) != 1:
            raise ValueError("Length of value for 'decchar' must be 1")

        c_detectorname = str_conversion(detectorname)
        c_decchar = ord(str_conversion(decchar))
        c_desc = str_conversion(desc)
        PY_SCIP_CALL( DECincludeDetector(
            self._scip, c_detectorname, c_decchar, c_desc, freqcallround, maxcallround, mincallround,
            freqcallroundoriginal, maxcallroundoriginal, mincallroundoriginal, priority, enabled, enabledfinishing,
            enabledpostprocessing, skip, usefulrecall, <DEC_DETECTORDATA*>detector, PyDetectorFree, PyDetectorInit,
            PyDetectorExit, PyDetectorPropagatePartialdec, PyDetectorFinishPartialdec, PyDetectorPostprocessPartialdec,
            PyDetectorSetParamAggressive, PyDetectorSetParamDefault, PyDetectorSetParamFast) )

        detector.model = <Model>weakref.proxy(self)
        detector.detectorname = detectorname
        Py_INCREF(detector)

    def listDetectors(self):
        """!@brief Lists all detectors that are currently included
        @return A list of strings of the detector names

        The detectors can be enabled or disabled using the appropriate methods by passing the name.

        @see setDetectorEnabled()
        @see setDetectorFinishingEnabled()
        @see setDetectorPostprocessingEnabled()
        """
        cdef int n_detectors = GCGconshdlrDecompGetNDetectors(self._scip)
        cdef DEC_DETECTOR** detectors = GCGconshdlrDecompGetDetectors(self._scip)

        return [DECdetectorGetName(detectors[i]).decode('utf-8') for i in range(n_detectors)]

    def setDetectorEnabled(self, detector_name, is_enabled=True):
        """!@brief Enables or disables a detector for detecting partial decompositions.
        @param detector_name The name of the detector.
        @param is_enabled Decides weather the detector should be enabled or diabled.

        This is a convenience method to access the boolean parameter "detection/detectors/<name>/enabled".

        @note Disabling a detector using this method is not enough to ensure that it will not run. In addition setDetectorFinishingEnabled() and setDetectorPostProcessingEnabled() have to be used.

        Use listDetectors() to obtain a list of all detectors.
        """
        # TODO test if detector_name exists
        self.setBoolParam("detection/detectors/{}/enabled".format(detector_name), is_enabled)

    def setDetectorFinishingEnabled(self, detector_name, is_enabled=True):
        """!@brief Enables or disables a detector for finishing partial decompositions.
        @param detector_name The name of the detector.
        @param is_enabled Decides weather the detector should be enabled or diabled.

        This is a convenience method to access the boolean parameter "detection/detectors/<name>/finishingenabled".

        @see setDetectorEnabled()
        """
        # TODO test if detector_name exists
        self.setBoolParam("detection/detectors/{}/finishingenabled".format(detector_name), is_enabled)

    def setDetectorPostprocessingEnabled(self, detector_name, is_enabled=True):
        """!@brief Enables or disables a detector for postprocessing partial decompositions.
        @param detector_name The name of the detector.
        @param is_enabled Decides weather the detector should be enabled or diabled.

        This is a convenience method to access the boolean parameter "detection/detectors/<name>/postprocessingenabled".

        @see setDetectorEnabled()
        """
        # TODO test if detector_name exists
        self.setBoolParam("detection/detectors/{}/postprocessingenabled".format(detector_name), is_enabled)

    def getMasterProb(self):
        """!@brief Provides access to the GCG master problem.
        @return An instance of scip#Model that represents the master problem.
        """
        cdef SCIP * master_prob = GCGgetMasterprob(self._scip)
        return Model.create(master_prob)

    def setGCGSeparating(self, setting):
        """!@brief Sets parameter settings of all separators
        @param setting: the parameter settings (SCIP_PARAMSETTING)
        """
        # GCG API is inconsistant with SCIP, SCIPsetSeparating
        PY_SCIP_CALL(GCGsetSeparators(self._scip, setting))

    def writeAllDecomps(self, directory="alldecompositions/", extension="dec", bool original=True, bool presolved=True, createDirectory=True):
        """!@brief Writes all decompositions to disk
        @param directory: A path to a folder where to store the decomposition files
        @param extension: Extension without a dot. Decides the output format. Use "dec" to output decomposition files
        @param createDirectory: Automatically create the directory specified in @p directory if it does not exist
        """
        if createDirectory:
            Path(directory).mkdir(exist_ok=True, parents=True)
        c_directory = str_conversion(directory)
        c_extension = str_conversion(extension)
        PY_SCIP_CALL( DECwriteAllDecomps(self._scip, c_directory, c_extension, original, presolved) )


cdef class PartialDecomposition:
    """!@brief class to manage partial decompositions

    each partialdec corresponds to one detprobdata which contains the problem information,
    there is one detprobdata for the original and the transformed problem.
    """
    cdef PARTIALDECOMP * thisptr
    cdef bool delete_thisptr

    cdef public dict _visualizations

    # Stores the objects used by the user for referencing the decomposition blocks.
    # The index in the list corresponds to GCG's internal block id number.
    cdef list py_block_id_map

    def __cinit__(self):
        self.thisptr = NULL
        self.delete_thisptr = True

    def __dealloc__(self):
        if self.delete_thisptr and self.thisptr != NULL:
            del self.thisptr

    @staticmethod
    cdef create(PARTIALDECOMP* thisptr):
        if thisptr == NULL:
            raise Warning("cannot create PartialDecomposition with PARTIALDECOMP* == NULL")
        new_PartialDecomposition = PartialDecomposition()
        new_PartialDecomposition.thisptr = thisptr
        new_PartialDecomposition._visualizations = {}
        new_PartialDecomposition.py_block_id_map = []
        new_PartialDecomposition.delete_thisptr = False
        return new_PartialDecomposition

    # def __init__(PartialDecomposition self, SCIP * scip, bool originalProblem):
    #     """@brief Standard constructor, creates empty partialdec with unique id
    #     @note initially, all conss and vars are open.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def __init__(PartialDecomposition self, PartialDecomposition partialdecToCopy):
    #     """!@brief copy constructor.
    #     """
    #     cdef PARTIALDECOMP * cpp_partialdecToCopy = partialdecToCopy.thisptr
    #     self.thisptr = new PARTIALDECOMP(cpp_partialdecToCopy)

    def __copy__(self):
        return PartialDecomposition.create(new PARTIALDECOMP(self.thisptr))

    def copy(self):
        return copy(self)

    def fixConsToMaster(PartialDecomposition self, Constraint cons):
        """!@brief fixes a Constraint to the master constraints
        @param cons scip#Constraint to add
        """
        self.thisptr.fixConsToMaster(cons.scip_cons)

    def fixConssToMaster(self, conss):
        """!@brief Fixes all constraints to the master constraints
        @param conss An iterable of scip#Constraint objects
        @see fixConsToMaster()
        """
        if not isinstance(conss, Iterable):
            raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(conss)))
        for cons in conss:
            self.fixConsToMaster(<Constraint?>cons)

    def fixConsToBlockId(PartialDecomposition self, Constraint cons, int block_id):
        """!@brief adds a constraint to a block
        @param cons scip#Constraint to add
        @param block_id id of block to add.
        @note The passed @p block_id has to by of type integer and use the internal numbering of blocks. Before
        calling this method, one has to ensure that the specified block exists. A more convenient alternative is
        fixConsToBlock().

        @see addBlock()
        """
        self.thisptr.fixConsToBlock(cons.scip_cons, block_id)

    def fixConssToBlockId(self, conss, int block_id):
        """!@brief adds all constraints to a block
        @param conss An iterable of scip#Constraint objects
        @param block_id id of block to add.
        @see fixConsToBlockId()
        @see fixConssToBlock()
        """
        if not isinstance(conss, Iterable):
            raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(conss)))
        for cons in conss:
            self.fixConsToBlockId(<Constraint?>cons, block_id)

    def fixConsToBlock(self, Constraint cons, object block):
        """!@brief adds a constraint to a block
        @param cons scip#Constraint to add
        @param block identifier of block to add. Can be any hashable Python object.
        @return the internal block_id assigned to the block.

        This method will automatically manage the blocks of the decomposition and create blocks if neccessary.
        The passed @p block identifiers will be mapped to internal block ids.

        To address the internal blocks directly, use @ref fixConsToBlockId().

        @see fixConssToBlock()
        """
        if block not in self.py_block_id_map:
            self.py_block_id_map.append(block)

        # Extend the number of blocks if we do not have enough. Do not automatically shrink the block count to avoid
        # interferance in case the user previously set a block count manually.
        if len(self.py_block_id_map) > self.getNBlocks():
            for i in range(len(self.py_block_id_map) - self.getNBlocks()):
                self.addBlock()

        cdef int block_id = self.py_block_id_map.index(block)
        self.fixConsToBlockId(cons, block_id)
        return block_id

    def fixConssToBlock(self, conss, object block):
        """!@brief adds all constraints to a block
        @param conss An iterable of scip#Constraint objects
        @param block identifier of block to add. Can be any hashable Python object.
        @return the internal block_id assigned to the block.

        @see fixConsToBlock()
        """
        if not isinstance(conss, Iterable):
            raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(conss)))
        cdef int block_id
        for cons in conss:
            # block_id is the same for all blocks
            block_id = self.fixConsToBlock(<Constraint?>cons, block)
        return block_id

    def getOpenconss(PartialDecomposition self):
        """!@brief Gets a vector containing constraint ids not assigned yet as vector
        @return returns a vector containing constraint ids not assigned yet as vector.
        """
        cdef vector[int] result = self.thisptr.getOpenconssVec()
        cdef DETPROBDATA * det_prob_data = self.thisptr.getDetprobdata()
        return [Constraint.create(det_prob_data.getCons(consIndex)) for consIndex in result]

    def setUsergiven(self, USERGIVEN value=COMPLETED_CONSTOMASTER):
        self.thisptr.setUsergiven(value)

    @property
    def classic_score(self):
        return self.thisptr.getClassicScore()

    @property
    def border_area_score(self):
        return self.thisptr.getBorderAreaScore()

    @property
    def max_white_score(self):
        return self.thisptr.getMaxWhiteScore()

    @property
    def max_for_white_score(self):
        return self.thisptr.getMaxForWhiteScore()

    @property
    def set_part_for_white_score(self):
        return self.thisptr.getSetPartForWhiteScore()

    @property
    def max_for_white_agg_score(self):
        return self.thisptr.getMaxForWhiteAggScore()

    @property
    def set_part_for_white_agg_score(self):
        return self.thisptr.getSetPartForWhiteAggScore()

    @property
    def benders_score(self):
        return self.thisptr.getBendersScore()

    @property
    def strong_decomp_score(self):
        return self.thisptr.getStrongDecompScore()

    # BEGIN AUTOGENERATED BLOCK

    def addBlock(PartialDecomposition self):
        """!@brief adds a block
        @returns the number (id) of the new block
        .
        """
        cdef int result = self.thisptr.addBlock()
        return result


    def addClockTime(PartialDecomposition self, double clocktime):
        """!@brief adds detection time of one detector

        incorporates the needed time of some detector in the detector chain.
        """
        cdef double cpp_clocktime = clocktime
        self.thisptr.addClockTime(cpp_clocktime)

    def addDecChangesFromAncestor(PartialDecomposition self, PartialDecomposition ancestor):
        """!@brief adds the statistical differences to an ancestor

        incorporates the changes from ancestor partialdec into the statistical data structures.
        """
        cdef PARTIALDECOMP * cpp_ancestor = ancestor.thisptr
        self.thisptr.addDecChangesFromAncestor(cpp_ancestor)

    def addDetectorChainInfo(PartialDecomposition self, decinfo):
        """!@brief add information about the detector chain

        adds a detectorchain information string to the corresponding vector
        (that carries information for each detector call)
        .
        """
        c_decinfo = str_conversion(decinfo)
        self.thisptr.addDetectorChainInfo(c_decinfo)

    def addNNewBlocks(PartialDecomposition self, int nnewblocks):
        """!@brief adds how many new blocks were introduced

        bookkeeping information: adds number of new blocks created by a detector added to detector chain.
        """
        cdef int cpp_nnewblocks = nnewblocks
        self.thisptr.addNNewBlocks(cpp_nnewblocks)

    def addPctConssFromFree(PartialDecomposition self, double pct):
        """!@brief adds percentage of closed constraints

        bookkeeping information: fraction of constraints that are not longer open for a detector added to detector chain.
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssFromFree(cpp_pct)

    def addPctConssToBlock(PartialDecomposition self, double pct):
        """!@brief adds percentage of constraints assigned to blocks

        bookkeeping information: adds fraction of constraints assigned to a block for a detector added to detector chain
        .
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssToBlock(cpp_pct)

    def addPctConssToBorder(PartialDecomposition self, double pct):
        """!@brief adds percentage of constraints assigned to border

        bookkeeping information: adds fraction of constraints assigned to the border for a detector added to detector chain.
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssToBorder(cpp_pct)

    def addPctVarsFromFree(PartialDecomposition self, double pct):
        """!@brief adds percentage of closed variables

        bookkeeping information: adds fraction of variables that are not longer open for a detector added to detector chain.
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsFromFree(cpp_pct)

    def addPctVarsToBlock(PartialDecomposition self, double pct):
        """!@brief adds percentage of variables assigned to blocks

        bookkeeping information: adds fraction of variables assigned to a block for a detector added to detector chain
        .
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsToBlock(cpp_pct)

    def addPctVarsToBorder(PartialDecomposition self, double pct):
        """!@brief adds percentage of variables assigned to border

        bookkeeping information: adds fraction of variables assigned to the border for a detector added to detector chain.
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsToBorder(cpp_pct)

    def alreadyAssignedConssToBlocks(PartialDecomposition self):
        """!@brief method to check if at least one constraint is assigned to some block
        @returns true iff at least one constraint is assigned to a block
        .
        """
        cdef bool result = self.thisptr.alreadyAssignedConssToBlocks()
        return result


    # def assignBorderFromConstoblock(PartialDecomposition self, SCIP_HASHMAP * constoblock, int givenNBlocks):
    #     """@brief assigns open conss to master

    #     assigns open constraints to master according to the cons assignment information given in constoblock hashmap
    #     @returns scip return code
    #     @note for conss assigned to blocks according to constoblock there is no assignment \see assignPartialdecFromConstoblock
    #     @note master assignment is indicated by assigning cons to index additionalNBlocks
    #     .
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def assignCurrentStairlinking(PartialDecomposition self):
        """!@brief assigns open vars to stairlinking if appropriate

        assigns open vars to stairlinking if they can be found in exactly two consecutive blocks
        @returns true iff at least one stairlinkingvar was assigned.
        """
        cdef bool result = self.thisptr.assignCurrentStairlinking()
        return result


    def assignOpenConssToMaster(PartialDecomposition self):
        """!@brief assigns open conss to master.
        """
        self.thisptr.assignOpenConssToMaster()

    # def assignPartialdecFromConstoblock(PartialDecomposition self, SCIP_HASHMAP * constoblock, int additionalNBlocks):
    #     """@brief assigns conss structure according to given hashmap

    #     adds blocks and assigns open conss to a new block or to master
    #     according to the cons assignment information given in constoblock hashmap
    #     @returns scip return code
    #     \see assignPartialdecFromConstoblockVector()
    #     @note master assignment is indicated by assigning cons to index additionalNBlocks
    #     .
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def assignPartialdecFromConstoblockVector(PartialDecomposition self, object constoblock, int additionalNBlocks):
        """!/*!
        @brief assigns conss structure according to given vector

        adds blocks and assigns open conss to a new block or to master
        according to the cons assignment information given in constoblock vector
        @returns scip return code
        \see  assignPartialdecFromConstoblock()
        @note master is indicated by assigning cons to index additionalNBlocks.
        """
        cdef vector[int] cpp_constoblock = constoblock
        cdef int cpp_additionalNBlocks = additionalNBlocks
        # TODO implement function
        raise NotImplementedError()

    def assignSmallestComponentsButOneConssAdjacency(PartialDecomposition self):
        """!@brief computes components by connectedness of conss and vars

        computes components corresponding to connectedness of conss and vars
        and assigns them accordingly (all but one of largest components)

        strategy: assigns all conss same block if they are connected
        two constraints are adjacent if there is a common variable

        @note this relies on the consadjacency structure of the detprobdata
        hence it cannot be applied in presence of linking variables.
        """
        self.thisptr.assignSmallestComponentsButOneConssAdjacency()

    def calcStairlinkingVars(PartialDecomposition self):
        """!@brief reassigns linking vars to stairlinkingvars if possible

        potentially reorders blocks for making a maximum number of linking vars stairlinking
        if all vars that connect exactly two blocks have a staircase structure, all of them become stairlinkingvars
        otherwise, the stairlinking assignment is done greedily
        @note precondition: partialdec does not have any stairlinking vars.
        """
        self.thisptr.calcStairlinkingVars()

    def checkAllConssAssigned(PartialDecomposition self):
        """!@brief checks if all conss are assigned

        returns true iff all constraints are assigned and deletes the vector open conss if so
        @return true iff all constraints are assigned
        .
        """
        cdef bool result = self.thisptr.checkAllConssAssigned()
        return result


    def checkConsistency(PartialDecomposition self):
        """!@brief Checks whether the assignments in the partialdec are consistent

        The following checks are performed:
        - check if nblocks is set appropriately
        - check for empty (row- and col-wise) blocks
        - every variable is assigned at most once
        - check if all not assigned variables are open vars
        - check if all open vars are not assigned
        - every constraint is assigned at most once
        - check if all not assigned constraints are open cons
        - check if all open conss are not assigned
        - check if the data structures are sorted
        - check if variables hitting a cons are either in the cons's block or border or still open
        @returns true iff the partialdec seems to be consistent
        .
        """
        cdef bool result = self.thisptr.checkConsistency()
        return result


    def complete(PartialDecomposition self):
        """!@brief assigns all open constraints and open variables trivially

        strategy: assigns all open conss and vars to blocks if they can be refined there, otherwise to the master

        @note partialdecomps should usually be completed by a detector, only use this function if you know what you are doing.
        """
        self.thisptr.complete()

    def completeByConnected(PartialDecomposition self):
        """!@brief assigns all open constraints and open variables

        strategy: assigns all conss and vars to the same block if they are connected,
        a cons and a var are adjacent if the var appears in the cons.
        """
        self.thisptr.completeByConnected()

    def completeByConnectedConssAdjacency(PartialDecomposition self):
        """!@brief assigns all open constraints and open variables

        strategy: assigns all conss and vars to the same block if they are connected
        a cons and a var are adjacent if the var appears in the cons
        \note this relies on the consadjacency structure of the detprobdata
        hence it cannot be applied in presence of linking variables.
        """
        self.thisptr.completeByConnectedConssAdjacency()

    def completeGreedily(PartialDecomposition self):
        """!@brief assigns all open constraints and open variables

        strategy: assigns a cons (and related vars) to a new block if possible,
        if not to an existing block if possible (by means of prior var assignments)
        and finally to master, if there does not exist such a block.
        """
        self.thisptr.completeGreedily()

    def removeMastercons(PartialDecomposition self, int consid):
        """!@brief removes the given cons from master.
        """
        cdef int cpp_consid = consid
        self.thisptr.removeMastercons(cpp_consid)

    def considerImplicits(PartialDecomposition self):
        """!@brief: assigns every open cons/var

        Assignments happen as follows:
        - to the respective block if it hits exactly one blockvar/blockcons and no open vars/conss
        - to master/linking if it hits blockvars/blockconss assigned to different blocks
        - and every cons to master that hits a master var
        - and every var to master if it does not hit any blockcons and has no open cons
        - leave the cons/variableopen if nothing from the above holds
        .
        """
        self.thisptr.considerImplicits()

    def copyPartitionStatistics(PartialDecomposition self, PartialDecomposition otherpartialdec):
        """!@brief copies the given partialdec's partition statistics

        @param otherpartialdec partialdec whose partition statistics are to be copied.
        """
        cdef PARTIALDECOMP * cpp_otherpartialdec = otherpartialdec.thisptr
        self.thisptr.copyPartitionStatistics(cpp_otherpartialdec)

    def deleteEmptyBlocks(PartialDecomposition self, bool variables):
        """!@brief deletes empty blocks and sets nblocks accordingly

        A block is considered to be empty if no constraint is assigned to it,
        variables in blocks with no constraints become open

        @param variables if true, then blocks with no constraints but at least one variable are considered to be nonempty.
        """
        cdef bool cpp_variables = variables
        self.thisptr.deleteEmptyBlocks(cpp_variables)

    def deleteOpencons(PartialDecomposition self, int opencons):
        """!@brief deletes a cons from list of open conss

        @param opencons id of the cons that is not considered open anymore.
        """
        cdef int cpp_opencons = opencons
        self.thisptr.deleteOpencons(cpp_opencons)

    def deleteOpenvar(PartialDecomposition self, int openvar):
        """!@brief deletes a var from the list of open vars

        @param openvar id of the var that is not considered open anymore.
        """
        cdef int cpp_openvar = openvar
        self.thisptr.deleteOpenvar(cpp_openvar)

    def displayInfo(PartialDecomposition self, int detailLevel):
        """!@brief displays the relevant information of the partialdec

        @param detailLevel pass a value that indicates how detailed the output should be:
        0: brief overview
        1: block and detector info
        2: cons and var assignments.
        """
        cdef int cpp_detailLevel = detailLevel
        self.thisptr.displayInfo(cpp_detailLevel)

    # def filloutBorderFromConstoblock(PartialDecomposition self, SCIP_HASHMAP * constoblock, int givenNBlocks):
    #     """@brief every constraint is either assigned to master or open

    #     Assignment happens according to the cons assignment information given in constoblock hashmap,
    #     variables are set accordingly
    #     @note precondition: no constraint or variable is already assigned to a block
    #     @return scip return code.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def filloutPartialdecFromConstoblock(PartialDecomposition self, SCIP_HASHMAP * constoblock, int givenNBlocks):
    #     """@brief assigns all conss to master or a block

    #     Assignment happens according to the cons assignment information given in constoblock hashmap

    #     @return scip return code
    #     calculates implicit variable assignment through cons assignment
    #     @note precondition: no cons or var is already assigned to a block and constoblock contains information for every cons.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def findVarsLinkingToMaster(PartialDecomposition self):
        """!@brief reassigns linking variables to master if appropriate

        Variables are reassigned as master if the variable only hits master conss.
        """
        self.thisptr.findVarsLinkingToMaster()

    def findVarsLinkingToStairlinking(PartialDecomposition self):
        """!@brief reassigns variables classified as linking to stairlinking if appropriate

        Variables are reassigned as master if the variable hits conss in exactly two consecutive
        blocks.
        """
        self.thisptr.findVarsLinkingToStairlinking()

    def getAncestorID(PartialDecomposition self, int ancestorindex):
        """!@brief gets partialdec id of given ancestor id
        @return partialdec id of given ancestor id.
        """
        cdef int cpp_ancestorindex = ancestorindex
        cdef int result = self.thisptr.getAncestorID(cpp_ancestorindex)
        return result


    def getAncestorList(PartialDecomposition self):
        """!@brief get ancestor ids as vector
        @return vector of ids of all ancestors id.
        """
        cdef vector[int] result = self.thisptr.getAncestorList()
        return result


    def setAncestorList(PartialDecomposition self, object newlist):
        """!set ancestor list directly
        @param newlist new list of ancestor ids.
        """
        cdef vector[int] cpp_newlist = newlist
        self.thisptr.setAncestorList(cpp_newlist)

    def removeAncestorID(PartialDecomposition self, int ancestorid):
        """!removes ancestor id from list.
        """
        cdef int cpp_ancestorid = ancestorid
        self.thisptr.removeAncestorID(cpp_ancestorid)

    def addAncestorID(PartialDecomposition self, int ancestor):
        """!adds ancestor id to back of list
        @param ancestor id of ancestor that is to be added.
        """
        cdef int cpp_ancestor = ancestor
        self.thisptr.addAncestorID(cpp_ancestor)

    def getBlocksForRep(PartialDecomposition self, int repid):
        """!@brief get a vector of block ids that are identical to block with id repid
        @param repid id of the representative block
        @return vector of block ids that are identical to block with id repid.
        """
        cdef int cpp_repid = repid
        cdef vector[int] result = self.thisptr.getBlocksForRep(cpp_repid)
        return result


    def getDetectorClockTime(PartialDecomposition self, int detectorchainindex):
        """!@brief returns the time that the detector related to the given detectorchainindex needed for detecting
        @return the clock time for the corresponding detector in the chain.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getDetectorClockTime(cpp_detectorchainindex)
        return result


    def getDetectorClockTimes(PartialDecomposition self):
        """!@brief returns a vector of the clock times that each detector needed that was involved in this partialdec
        @return vector of the clock times.
        """
        cdef vector[double] result = self.thisptr.getDetectorClockTimes()
        return result


    def getBlockConss(PartialDecomposition self, int block):
        """!@brief returns array containing constraints assigned to a block
        @param block id of the block the constraint indices are returned
        @return array containing constraints assigned to a block.
        @note This calls the corresponding `getConssForBlock()` method. For consistancy, the method was renamed in the Python API.
        """
        cdef int cpp_block = block
        cdef vector[int] result = self.thisptr.getConssForBlock(cpp_block)
        return [self.getDetprobdata().getCons(c) for c in result]


    # def getDetectorchain(PartialDecomposition self):
    #     """@brief returns detector chain as vector of detector pointers
    #     @return detector chain as array of detector pointers.
    #     """
    #     cdef vector[DEC_DETECTOR *] result = self.thisptr.getDetectorchain()
    #     return result


    def getFinishedByFinisher(PartialDecomposition self):
        """!@brief returns true iff this partialdec was finished by finishPartialdec() method of a detector
        @return true iff this partialdec was finished by finishPartialdec() method of a detector.
        """
        cdef bool result = self.thisptr.getFinishedByFinisher()
        return result


    def getHashValue(PartialDecomposition self):
        """!@brief returns the calculated hash value of this partialdec
        @return the calculated hash value of this partialdec.
        """
        cdef unsigned long result = self.thisptr.getHashValue()
        return result


    def getID(PartialDecomposition self):
        """!@brief returns the unique id of the partialdec
        @return the unique id of the partialdec.
        """
        cdef int result = self.thisptr.getID()
        return result


    def getLinkingvars(PartialDecomposition self):
        """!@brief returns array containing all linking vars indices
        @return vector containing all linking vars indices
        @note when accessed it is supposed to be sorted.
        """
        cdef vector[int] result = self.thisptr.getLinkingvars()
        return result


    def getMasterconss(PartialDecomposition self):
        """!@brief Gets array containing all master conss indices
        @return array containing all master conss indices
        @note when accessed it is supposed to be sorted.
        """
        cdef vector[int] result = self.thisptr.getMasterconss()
        return [self.getDetprobdata().getCons(c) for c in result]


    def getMastervars(PartialDecomposition self):
        """!@brief Gets array containing all master vars indices

        master vars hit only constraints in the master, aka static variables
        @return array containing all master vars indices.
        """
        cdef vector[int] result = self.thisptr.getMastervars()
        return result


    def getNCoeffsForBlock(PartialDecomposition self, int blockid):
        """!@brief Gets the number of nonzero coeffs in a certain block
        @param blockid of the block the number of nozerors are requested for
        @return number of nonzero coeffs in a certain block.
        """
        cdef int cpp_blockid = blockid
        cdef int result = self.thisptr.getNCoeffsForBlock(cpp_blockid)
        return result


    def getNCoeffsForMaster(PartialDecomposition self):
        """!Gets the number of nonzero coeffs in master
        @return the number of nonzero coeffs in master.
        """
        cdef int result = self.thisptr.getNCoeffsForMaster()
        return result


    # def getScore(PartialDecomposition self, SCORETYPE type):
    #     """@brief returns the score of the partialdec (depending on used scoretype)
    #     @param type the scoretype
    #     @return the score
    #     @see enum scoretype in cons_decomp.

    #     h
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def hasSetppccardMaster(PartialDecomposition self):
        """!@brief checks if all master constraints set partitioning, set packing, set cover, or cardinality constraints
        @return TRUE iff all master constraints set partitioning, set packing, set cover, or cardinality constraints.
        """
        cdef unsigned int result = self.thisptr.hasSetppccardMaster()
        return result


    def hasSetppcMaster(PartialDecomposition self):
        """!@brief checks iff all master constraints set partitioning, set packing, or set cover constraints
        @return TRUE iff all master constraints set partitioning, set packing, or set cover.
        """
        cdef unsigned int result = self.thisptr.hasSetppcMaster()
        return result


    def hasSetppMaster(PartialDecomposition self):
        """!@brief checks iff all master constraints set partitioning, or set packing constraints
        @return TRUE iff all master constraints set partitioning, or set packing constraints.
        """
        cdef unsigned int result = self.thisptr.hasSetppMaster()
        return result


    def getUsergiven(PartialDecomposition self):
        """!@brief Gets the USERGIVEN status of this partialdecs
        @return the USERGIVEN status of this partialdecs
        @see enum USERGIVEN.
        """
        # TODO implement function
        raise NotImplementedError()

    def getNAncestors(PartialDecomposition self):
        """!@brief Gets number of ancestor partialdecs
        @return number of ancestor partialdecs.
        """
        cdef int result = self.thisptr.getNAncestors()
        return result


    def getNBlocks(PartialDecomposition self):
        """!@brief Gets the number of blocks
        @return number of blocks.
        """
        cdef int result = self.thisptr.getNBlocks()
        return result


    def getNConss(PartialDecomposition self):
        """!@brief Gets the number of constraints
        @return number of constraints.
        """
        cdef int result = self.thisptr.getNConss()
        return result


    def getNConssForBlock(PartialDecomposition self, int block):
        """!@brief Gets size of the vector containing conss assigned to a block
        @param block id of the block the number of constraints is asked for
        @return size of the vector containing conss assigned to a block.
        """
        cdef int cpp_block = block
        cdef int result = self.thisptr.getNConssForBlock(cpp_block)
        return result


    def getDetectorchainInfo(PartialDecomposition self):
        """!@brief Gets the detectorchain info vector
        @return detectorchain info vector.
        """
        cdef vector[string] result = self.thisptr.getDetectorchainInfo()
        return result


    def getNDetectors(PartialDecomposition self):
        """!@brief Gets the number of detectors the partialdec is propagated by
        @return number of detectors the partialdec is propagated by.
        """
        cdef int result = self.thisptr.getNDetectors()
        return result


    def getNLinkingvars(PartialDecomposition self):
        """!@brief Gets size of the vector containing linking vars
        @return size of the vector containing linking vars.
        """
        cdef int result = self.thisptr.getNLinkingvars()
        return result


    def getNMasterconss(PartialDecomposition self):
        """!@brief Gets size of the vector containing master conss
        @returns size of the vector containing master conss.
        """
        cdef int result = self.thisptr.getNMasterconss()
        return result


    def getNMastervars(PartialDecomposition self):
        """!@brief Gets size of the vector containing master vars

        master vars hit only constraints in the master
        @return size of the vector containing master vars.
        """
        cdef int result = self.thisptr.getNMastervars()
        return result


    def getNNewBlocks(PartialDecomposition self, int detectorchainindex):
        """!@brief Gets number of blocks a detector added

        @return number of blocks a detector added.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef int result = self.thisptr.getNNewBlocks(cpp_detectorchainindex)
        return result


    def getNNewBlocksVector(PartialDecomposition self):
        """!@brief gets number of blocks the detectors in the detectorchain added
        @return number of blocks the detectors in the detectorchain added.
        """
        cdef vector[int] result = self.thisptr.getNNewBlocksVector()
        return result


    def getNTotalStairlinkingvars(PartialDecomposition self):
        """!@brief Gets total number of stairlinking vars
        @return total number of stairlinking vars.
        """
        cdef int result = self.thisptr.getNTotalStairlinkingvars()
        return result


    def getNOpenconss(PartialDecomposition self):
        """!@brief Gets size of vector containing constraints not assigned yet
        @return returns size of vector containing constraints not assigned yet.
        """
        cdef int result = self.thisptr.getNOpenconss()
        return result


    def getNOpenvars(PartialDecomposition self):
        """!@brief Gets size of vector containing variables not assigned yet
        @return size of vector containing variables not assigned yet.
        """
        cdef int result = self.thisptr.getNOpenvars()
        return result


    def getNReps(PartialDecomposition self):
        """!@brief Gets the number of blockrepresentatives
        @return the number of blockrepresentatives.
        """
        cdef int result = self.thisptr.getNReps()
        return result


    def getNStairlinkingvars(PartialDecomposition self, int block):
        """!@brief Gets size of the vector containing stairlinking vars
        @param block id of the block the size of the stairlinking vector is asked for
        @return size of the vector containing stairlinking vars.
        """
        cdef int cpp_block = block
        cdef int result = self.thisptr.getNStairlinkingvars(cpp_block)
        return result


    def getNVars(PartialDecomposition self):
        """!@brief Gets number of vars
        @return number of vars.
        """
        cdef int result = self.thisptr.getNVars()
        return result


    def getNVarsForBlock(PartialDecomposition self, int block):
        """!@brief Gets size of the vector containing vars assigned to a block
        @param block id of the block the number of variables is asked for
        @return size of the vector containing vars assigned to a block.
        """
        cdef int cpp_block = block
        cdef int result = self.thisptr.getNVarsForBlock(cpp_block)
        return result


    def getNVarsForBlocks(PartialDecomposition self):
        """!@brief Gets overall number of vars assigned to a block
        @return number of vars that are assigned to any block.
        """
        cdef int result = self.thisptr.getNVarsForBlocks()
        return result


    # def getOpenconss(PartialDecomposition self):
    #     """@brief Gets array containing constraints not assigned yet
    #     @return array containing constraints not assigned yet.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def getOpenconssVec(PartialDecomposition self):
        """!@brief Gets a vector containing constraint ids not assigned yet as vector
        @return returns a vector containing constraint ids not assigned yet as vector.
        """
        cdef vector[int] result = self.thisptr.getOpenconssVec()
        return result


    def getOpenvars(PartialDecomposition self):
        """!@brief Gets array containing variables not assigned yet
        @return returns array containing variables not assigned yet.
        """
        # TODO implement function
        raise NotImplementedError()

    def getOpenvarsVec(PartialDecomposition self):
        """!Gets array containing variables not assigned yet as vector
        @return array containing variables not assigned yet as vector.
        """
        cdef vector[int] result = self.thisptr.getOpenvarsVec()
        return result


    def getPctVarsToBorder(PartialDecomposition self, int detectorchainindex):
        """!@brief Gets fraction of variables assigned to the border for a detector

        @return fraction of variables assigned to the border for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsToBorder(cpp_detectorchainindex)
        return result


    def getPctVarsToBorderVector(PartialDecomposition self):
        """!@brief Gets fraction of variables assigned to the border for detectors in detectorchain
        @return vector of fractions of variables assigned to the border for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctVarsToBorderVector()
        return result


    def getPctVarsToBlock(PartialDecomposition self, int detectorchainindex):
        """!@brief Gets fraction of variables assigned to a block for a detector

        @return fraction of variables assigned to a block for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsToBlock(cpp_detectorchainindex)
        return result


    def getPctVarsToBlockVector(PartialDecomposition self):
        """!@brief returns fraction of variables assigned to a block for detectors in detectorchain
        @return vector of fractions of variables assigned to a block for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctVarsToBlockVector()
        return result


    def getPctVarsFromFree(PartialDecomposition self, int detectorchainindex):
        """!@brief Gets fraction of variables that are not longer open for a detector

        @return index of the detector in the detectorchain.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsFromFree(cpp_detectorchainindex)
        return result


    def getPctVarsFromFreeVector(PartialDecomposition self):
        """!@brief Gets fraction of variables that are not longer open for detectors in detectorchain
        @return vector or fractions of variables that are not longer open for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctVarsFromFreeVector()
        return result


    def getPctConssToBorder(PartialDecomposition self, int detectorchainindex):
        """!@brief Gets fraction of constraints assigned to the border for a detector
        @return returns fraction of constraints assigned to the border for a detector
        /
        /**.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssToBorder(cpp_detectorchainindex)
        return result


    def getPctConssToBorderVector(PartialDecomposition self):
        """!@brief Gets fraction of constraints assigned to the border for detectors in detectorchain
        @return vector of fractions of constraints assigned to the border for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctConssToBorderVector()
        return result


    def getPctConssToBlock(PartialDecomposition self, int detectorchainindex):
        """!@brief Gets fraction of constraints assigned to a block for a detector
        @return fraction of constraints assigned to a block for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssToBlock(cpp_detectorchainindex)
        return result


    def getPctConssToBlockVector(PartialDecomposition self):
        """!@brief Gets fraction of constraints assigned to a block for detectors in detectorchain
        @return vector of fractions of constraints assigned to a block for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctConssToBlockVector()
        return result


    def getPctConssFromFree(PartialDecomposition self, int detectorchainindex):
        """!@brief Gets fraction of constraints that are not longer open for a detector
        @return fraction of constraints that are not longer open for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssFromFree(cpp_detectorchainindex)
        return result


    def getPctConssFromFreeVector(PartialDecomposition self):
        """!@brief Gets fraction of constraints that are not longer open for detectors in detectorchain
        @return vector of fractions of constraints that are not longer open for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctConssFromFreeVector()
        return result


    def getRepForBlock(PartialDecomposition self, int blockid):
        """!@brief Gets index of the representative block for a block, this might be blockid itself
        @param blockid id of the block the representative is asked for
        @return index of the representative block for a block, this might be blockid itself.
        """
        cdef int cpp_blockid = blockid
        cdef int result = self.thisptr.getRepForBlock(cpp_blockid)
        return result


    def getRepVarmap(PartialDecomposition self, int repid, int blockrepid):
        """!@brief Gets the represenation varmap

        Var map is vector for represenative repid and the blockrepid-th block that is represented by repid
        @param repid id of representative
        @param blockrepid id of block
        @return the represenation varmap as vector for represenative repid and the blockrepid-th block that is represented by repid.
        """
        cdef int cpp_repid = repid
        cdef int cpp_blockrepid = blockrepid
        cdef vector[int] result = self.thisptr.getRepVarmap(cpp_repid, cpp_blockrepid)
        return result


    def getDetprobdata(PartialDecomposition self):
        """!@brief Gets the corresponding detprobdata
        @return corresponding detprobdata.
        """
        cdef DETPROBDATA * result = self.thisptr.getDetprobdata()
        return DetProbData.create(result)

    def getStairlinkingvars(PartialDecomposition self, int block):
        """!@brief Gets array containing stairlinking vars,
        @note if a stairlinking variable links block i and i+1 it is only stored in vector of block i
        @param block id of the block the stairlinking variable varctor is asked for
        @return array containing stairlinking vars,.
        """
        cdef int cpp_block = block
        # TODO implement function
        raise NotImplementedError()

    def getVarsForBlock(PartialDecomposition self, int block):
        """!@brief Gets array containing vars of a block
        @param block id of the block the vars are requested for
        @return returns array containing vars of a block.
        """
        cdef int cpp_block = block
        cdef vector[int] result = self.thisptr.getVarsForBlock(cpp_block)
        return result


    def getVarProbindexForBlock(PartialDecomposition self, int varid, int block):
        """!@brief  Gets index in variables array of a block for a variable
        @param varid the id of the variable the index
        @param block the corresponding block id
        @return  returns index in variables array of a block for a variable.
        """
        cdef int cpp_varid = varid
        cdef int cpp_block = block
        cdef int result = self.thisptr.getVarProbindexForBlock(cpp_varid, cpp_block)
        return result


    def isComplete(PartialDecomposition self):
        """!@brief Gets whether this partialdec is complete,
        i.

        e. it has no more open constraints and variables
        @return TRUE iff this partialdec is complete
        """
        cdef bool result = self.thisptr.isComplete()
        return result


    def isConsMastercons(PartialDecomposition self, int cons):
        """!@brief Gets whether the cons is a master cons
        @param cons id of ccons to check if it is master constraint
        @return true iff the cons is a master cons.
        """
        cdef int cpp_cons = cons
        cdef bool result = self.thisptr.isConsMastercons(cpp_cons)
        return result


    def isConsOpencons(PartialDecomposition self, int cons):
        """!@brief Gets whether the cons is an open cons
        @param cons id of cons to check
        @return true iff the cons is an open cons.
        """
        cdef int cpp_cons = cons
        cdef bool result = self.thisptr.isConsOpencons(cpp_cons)
        return result


    def isAssignedToOrigProb(PartialDecomposition self):
        """!@brief Gets whether the partialdec is from the presolved problem
        @return true iff the partialdec is from the presolved problem.
        """
        cdef bool result = self.thisptr.isAssignedToOrigProb()
        return result


    @property
    def isSelected(PartialDecomposition self):
        """!Gets whether the partialdec is currently selected in explore menue
        @return true iff the partialdec is currently selected in explore menue.
        """
        cdef bool result = self.thisptr.isSelected()
        return result

    @isSelected.setter
    def isSelected(PartialDecomposition self, bool selected):
        """!@brief set the selection status of this partialdecs
        @param selected whether the partialdec is selected.
        """
        cdef bool cpp_selected = selected
        self.thisptr.setSelected(cpp_selected)

    # def isEqual(PartialDecomposition self, PartialDecomposition otherpartialdec, unsigned int * isequal, bool sortpartialdecs):
    #     """@brief method to check whether this partialdec is equal to a given other partialdec ( \see  isEqual(PARTIALDECOMP*))

    #     @return scip return code.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def isPropagatedBy(PartialDecomposition self, DEC_DETECTOR * detector):
    #     """@brief Gets whether this partialdec was propagated by specified detector
    #     @param detector pointer to detector to check for
    #     @return true iff this partialdec was propagated by detectorID.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def isTrivial(PartialDecomposition self):
        """!@brief Gets whether this partialdec is considered to be trivial

        PARTIALDECOMP is considered trivial if all conss are in one block, all conss are in border,
        all variables linking or mastervars, or all constraints and variables are open
        @return true iff this partialdec is considered to be trivial.
        """
        cdef bool result = self.thisptr.isTrivial()
        return result


    def isVarBlockvarOfBlock(PartialDecomposition self, int var, int block):
        """!@brief Checks whether the var is assigned to the block
        @param var id of var to check
        @param block id of block to check
        @return true iff the var is assigned to the block.
        """
        cdef int cpp_var = var
        cdef int cpp_block = block
        cdef bool result = self.thisptr.isVarBlockvarOfBlock(cpp_var, cpp_block)
        return result


    def isVarLinkingvar(PartialDecomposition self, int var):
        """!@brief Checks whether the var is a linking var
        @param var id of var to check
        @return true iff the var is a linking var.
        """
        cdef int cpp_var = var
        cdef bool result = self.thisptr.isVarLinkingvar(cpp_var)
        return result


    def isVarMastervar(PartialDecomposition self, int var):
        """!@brief Checks whether the var is a master var
        @param var id of var to check
        @return true iff the var is a master var.
        """
        cdef int cpp_var = var
        cdef bool result = self.thisptr.isVarMastervar(cpp_var)
        return result


    def isVarOpenvar(PartialDecomposition self, int var):
        """!@brief Checks whether the var is an open var
        @param var id of var to check
        @return true iff the var is an open var
        /
        /**.
        """
        cdef int cpp_var = var
        cdef bool result = self.thisptr.isVarOpenvar(cpp_var)
        return result


    def isVarStairlinkingvar(PartialDecomposition self, int var):
        """!@brief Checks whether the var is a stairlinking var
        @param var id of var to check
        @return true iff the var is a stairlinking var.
        """
        cdef int cpp_var = var
        cdef bool result = self.thisptr.isVarStairlinkingvar(cpp_var)
        return result


    def isVarStairlinkingvarOfBlock(PartialDecomposition self, int var, int block):
        """!@brief Checks whether the var is a stairlinkingvar of a specified block
        @param var id of var to check if it is a stairlinking variable hitting specified block
        @param block id of block to check
        @return true iff the var is a stairlinkingvar of a specified block.
        """
        cdef int cpp_var = var
        cdef int cpp_block = block
        cdef bool result = self.thisptr.isVarStairlinkingvarOfBlock(cpp_var, cpp_block)
        return result


    # def printPartitionInformation(PartialDecomposition self, SCIP * givenscip, FILE * file):
    #     """@brief prints partition information as described in \see cls reader
    #     @param givenscip scip data structure
    #     @param file output file.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def refineToBlocks(PartialDecomposition self):
        """!@brief refine partialdec with focus on blocks

        strategy: assigns open conss and vars if they can be found in blocks
        (without respect to open vars and conss  @see assignHittingOpenconss(), @see assignHittingOpenvars())
        @note partialdec might be not complete.
        """
        self.thisptr.refineToBlocks()

    def refineToMaster(PartialDecomposition self):
        """!@brief refine partialdec with focus on master

        strategy: do obvious ( @see considerImplicits()) assignments and
        assign other conss and vars to master if possible (@see assignOpenPartialHittingToMaster()).
        """
        self.thisptr.refineToMaster()

    def setConsPartitionStatistics(PartialDecomposition self, int detectorchainindex, ConsPart partition, object consclassesmaster):
        """!@brief registers statistics for a used conspartition.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef ConsPartition * cpp_partition = partition.thisptr
        cdef vector[int] cpp_consclassesmaster = consclassesmaster
        self.thisptr.setConsPartitionStatistics(cpp_detectorchainindex, cpp_partition, cpp_consclassesmaster)

    def setConsToBlock(PartialDecomposition self, int consToBlock, int block):
        """!@brief adds a constraint to a block, does not delete this cons from list of open conss
        @param consToBlock id of cons to add
        @param block id of block to add.
        """
        cdef int cpp_consToBlock = consToBlock
        cdef int cpp_block = block
        self.thisptr.setConsToBlock(cpp_consToBlock, cpp_block)

    # def fixConsToBlock(PartialDecomposition self, int cons, int block):
    #     """@brief adds a constraint to a block
    #     @param cons id of cons to add
    #     @param block id of block to add.
    #     """
    #     cdef int cpp_cons = cons
    #     cdef int cpp_block = block
    #     self.thisptr.fixConsToBlock(cpp_cons, cpp_block)

    def setConsToMaster(PartialDecomposition self, int consToMaster):
        """!@brief adds a constraint to the master constraints, does not delete this cons from list of open conss
        @param consToMaster id of cons to add.
        """
        cdef int cpp_consToMaster = consToMaster
        self.thisptr.setConsToMaster(cpp_consToMaster)

    # def fixConsToMaster(PartialDecomposition self, int cons):
    #     """@brief fixes a constraint to the master constraints
    #     @param cons id of cons to add
    #     @warning This method modifies the vector PARTIALDECOMP::openconss! Hence, any kind of iterator might be invalid afterwards!.
    #     """
    #     cdef int cpp_cons = cons
    #     self.thisptr.fixConsToMaster(cpp_cons)

    # def setDetectorchain(PartialDecomposition self, object givenDetectorChain):
    #     """@brief sets the detectorchain with the given vector of detector pointers
    #     @param givenDetectorChain vector of detector pointers.
    #     """
    #     cdef vector[DEC_DETECTOR *] cpp_givenDetectorChain = givenDetectorChain
    #     self.thisptr.setDetectorchain(cpp_givenDetectorChain)

    # def setDetectorPropagated(PartialDecomposition self, DEC_DETECTOR * detector):
    #     """@brief sets partialdec to be propagated by a detector
    #     @param detector pointer to detector that is registered for this partialdec.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def setDetectorFinished(PartialDecomposition self, DEC_DETECTOR * detector):
    #     """@brief sets detector that finished the partialdec
    #     @param detector pointer to detector that has finished this partialdecs.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def setDetectorFinishedOrig(PartialDecomposition self, DEC_DETECTOR * detectorID):
    #     """@brief sets detector that finished the partialdec in the original problem
    #     @param detectorID pointer to detector that has finished this partialdecs
    #     @note does not add the detector to the detectorchain and does not modify partition statistics.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def setFinishedByFinisher(PartialDecomposition self, bool finished):
        """!@brief sets whether this partialdec was finished by a finishing detector
        @param finished is this partialdecs finished by a finishing detector.
        """
        cdef bool cpp_finished = finished
        self.thisptr.setFinishedByFinisher(cpp_finished)

    def setFinishedByFinisherOrig(PartialDecomposition self, bool finished):
        """!@brief sets whether this partialdec was finished by a finishing detector in the original problem

        (in case this partialdec was translated)
        @param finished was this partialdecs finished by a finishing detector in orig.
        """
        cdef bool cpp_finished = finished
        self.thisptr.setFinishedByFinisherOrig(cpp_finished)

    def setNBlocks(PartialDecomposition self, int nblocks):
        """!@brief sets number of blocks, only increasing number allowed
        @param nblocks new number of blocks.
        """
        cdef int cpp_nblocks = nblocks
        self.thisptr.setNBlocks(cpp_nblocks)

    def setStemsFromOrig(PartialDecomposition self, bool fromorig):
        """!@brief sets whether this partialdec stems from an orig problem partialdec
        @param fromorig has this partialdec ancestors from the orig problem.
        """
        cdef bool cpp_fromorig = fromorig
        self.thisptr.setStemsFromOrig(cpp_fromorig)

    # def setUsergiven(PartialDecomposition self, cpp.USERGIVEN usergiven):
    #     """@brief sets whether this partialdec is user given
    #     @param usergiven is this partialdec user given.
    #     """
    #     self.thisptr.setUsergiven(usergiven)

    def setVarPartitionStatistics(PartialDecomposition self, int detectorchainindex, VarPart partition, object varclasseslinking, object varclassesmaster):
        """!@brief registers statistics for a used varpartition.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef VarPartition * cpp_partition = partition.thisptr
        cdef vector[int] cpp_varclasseslinking = varclasseslinking
        cdef vector[int] cpp_varclassesmaster = varclassesmaster
        self.thisptr.setVarPartitionStatistics(cpp_detectorchainindex, cpp_partition, cpp_varclasseslinking, cpp_varclassesmaster)

    def setVarToBlock(PartialDecomposition self, int varToBlock, int block):
        """!@brief adds a variable to the linking variables, does not delete this var from list of open vars
        @param varToBlock id of var to be added
        @param block id of block to be added.
        """
        cdef int cpp_varToBlock = varToBlock
        cdef int cpp_block = block
        self.thisptr.setVarToBlock(cpp_varToBlock, cpp_block)

    def fixVarToBlock(PartialDecomposition self, int var, int block):
        """!@brief adds a variable to the linking variables
        @param var id of var to be added
        @param block id of block to be added.
        """
        cdef int cpp_var = var
        cdef int cpp_block = block
        self.thisptr.fixVarToBlock(cpp_var, cpp_block)

    def setVarToLinking(PartialDecomposition self, int varToLinking):
        """!@brief adds a variable to the linking variables, does not delete this var from list of open vars
        @param varToLinking var to be set to linking.
        """
        cdef int cpp_varToLinking = varToLinking
        self.thisptr.setVarToLinking(cpp_varToLinking)

    def fixVarToLinking(PartialDecomposition self, int var):
        """!@brief adds a variable to the linking variables
        @param var var to be set to linking.
        """
        cdef int cpp_var = var
        self.thisptr.fixVarToLinking(cpp_var)

    def setVarToMaster(PartialDecomposition self, int varToMaster):
        """!@brief adds a variable to the master variables, does not delete this var from list of open vars

        master variables hit only constraints in the master.
        """
        cdef int cpp_varToMaster = varToMaster
        self.thisptr.setVarToMaster(cpp_varToMaster)

    def fixVarToMaster(PartialDecomposition self, int var):
        """!@brief adds a variable to the master variables

        master variables hit only constraints in the master.
        """
        cdef int cpp_var = var
        self.thisptr.fixVarToMaster(cpp_var)

    def setVarToStairlinking(PartialDecomposition self, int varToStairLinking, int block1, int block2):
        """!@brief adds a variable to the stairlinking variables, does not delete this var from list of open vars
        @param varToStairLinking id of variable to be added
        @param block1 id of block one
        @param block2 id of block two
        @note stairlinking variables are only registered in block with smaller index.
        """
        cdef int cpp_varToStairLinking = varToStairLinking
        cdef int cpp_block1 = block1
        cdef int cpp_block2 = block2
        self.thisptr.setVarToStairlinking(cpp_varToStairLinking, cpp_block1, cpp_block2)

    def fixVarToStairlinking(PartialDecomposition self, int var, int firstblock):
        """!@brief adds a variable to the stairlinking variables
        @param var id of variable to be added
        @param firstblock stairlinking variables hit exactly two consecutive blocks, this is the index of the first of these blocks
        @note stairlinking variables are only registered in block with smaller index.
        """
        cdef int cpp_var = var
        cdef int cpp_firstblock = firstblock
        self.thisptr.fixVarToStairlinking(cpp_var, cpp_firstblock)

    def fixConsToBlockByName(PartialDecomposition self, consname, int blockid):
        """!@brief assigns a constraint by name to a block
        @see fixConsToBlock
        @returns true iff successful.
        """
        c_consname = str_conversion(consname)
        cdef int cpp_blockid = blockid
        cdef bool result = self.thisptr.fixConsToBlockByName(c_consname, cpp_blockid)
        return result


    def fixVarToBlockByName(PartialDecomposition self, varname, int blockid):
        """!@brief assigns a variable by name to a block
        @see fixVarToBlock
        @returns true iff successful.
        """
        c_varname = str_conversion(varname)
        cdef int cpp_blockid = blockid
        cdef bool result = self.thisptr.fixVarToBlockByName(c_varname, cpp_blockid)
        return result


    def fixConsToMasterByName(PartialDecomposition self, consname):
        """!@brief assgins a constraint by name as master
        @see fixConsToMaster
        @returns true iff successful.
        """
        c_consname = str_conversion(consname)
        cdef bool result = self.thisptr.fixConsToMasterByName(c_consname)
        return result


    def fixVarToMasterByName(PartialDecomposition self, varname):
        """!@brief assigns a variable with given name as master
        @see fixVarToMaster
        @returns true iff successful.
        """
        c_varname = str_conversion(varname)
        cdef bool result = self.thisptr.fixVarToMasterByName(c_varname)
        return result


    def fixVarToLinkingByName(PartialDecomposition self, varname):
        """!@brief assigns a variable by name to the linking variables
        @see fixVarToLinking
        @returns true iff successful.
        """
        c_varname = str_conversion(varname)
        cdef bool result = self.thisptr.fixVarToLinkingByName(c_varname)
        return result


    def showVisualisation(PartialDecomposition self):
        """!@brief generates and opens a gp visualization of the partialdec
        @see visual/pdfreader and
        @note linux only.
        """
        self.thisptr.showVisualisation()

    # def generateVisualisation(PartialDecomposition self, filename, outname, GP_OUTPUT_FORMAT outputformat):
    #     """@brief generates a visualization of the partialdec using gnuplot
    #     @param filename Path where to store the gp file
    #     @param outname Path at which gnuplot will output its result
    #     @param outputformat The format of the gnuplot output.

    #     Should match the file extension of outname
    #     @note linux only, requires gnuplot
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def writeVisualisationFile(PartialDecomposition self, filename, outname, GP_OUTPUT_FORMAT outputformat):
    #     """@brief writes a gp visualization of the partialdec to a file
    #     @param filename Path where to store the gp file
    #     @param outname Path at which gnuplot will output its result
    #     @param outputformat The format of the gnuplot output.

    #     Should match the file extension of outname
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def shouldCompletedByConsToMaster(PartialDecomposition self):
        """!@brief Checks whether this partialdec is a userpartialdec that should be completed

        the completion should be done by setting unspecified constraints to master
        @return TRUE iff this partialdec is a userpartialdec that should be completed.
        """
        cdef unsigned int result = self.thisptr.shouldCompletedByConsToMaster()
        return result


    def sort(PartialDecomposition self):
        """!@brief sorts the vars and conss data structures by their indices
        @returns true if the internal order of variables or constraints changed.
        """
        cdef bool result = self.thisptr.sort()
        return result


    def setPctConssToBlockVector(PartialDecomposition self, object newvector):
        """!@brief set statistical vector of fractions of constraints set to blocks per involved detector
        @param newvector vector of fractions of constraints set to blocks per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssToBlockVector(cpp_newvector)

    def setPctConssFromFreeVector(PartialDecomposition self, object newvector):
        """!@brief set statistical vector of fractions of constraints that are not longer open per involved detector
        @param newvector vector of fractions of constraints that are not longer open per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssFromFreeVector(cpp_newvector)

    def setPctConssToBorderVector(PartialDecomposition self, object newvector):
        """!@brief set statistical vector of fractions of constraints assigned to the border per involved detector
        @param newvector vector of fractions of constraints assigned to the border per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssToBorderVector(cpp_newvector)

    def setPctVarsToBorderVector(PartialDecomposition self, object newvector):
        """!@brief set statistical vector of fraction of variables assigned to the border per involved detector
        @param newvector vector of fractions of variables assigned to the border per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsToBorderVector(cpp_newvector)

    def setPctVarsToBlockVector(PartialDecomposition self, object newvector):
        """!@brief set statistical vector of fractions of variables assigned to a block per involved detector
        @param newvector vector of fractions of variables assigned to a block per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsToBlockVector(cpp_newvector)

    def setPctVarsFromFreeVector(PartialDecomposition self, object newvector):
        """!@brief set statistical vector of variables that are not longer open per involved detector
        @param newvector vector of fractions of variables that are not longer open per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsFromFreeVector(cpp_newvector)

    def setDetectorClockTimes(PartialDecomposition self, object newvector):
        """!@brief set statistical vector of the times that the detectors needed for detecting per involved detector
        @param newvector vector of the times that the detectors needed for detecting per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setDetectorClockTimes(cpp_newvector)

    @property
    def classicScore(PartialDecomposition self):
        """!@brief gets the classic score

        @note -1 iff not calculated yet, \see GCGconshdlrDecompCalcClassicScore
        @returns border area score.
        """
        cdef double result = self.thisptr.getClassicScore()
        return result


    @classicScore.setter
    def classicScore(PartialDecomposition self, double score):
        """!@brief set the classic score.
        """
        cdef double cpp_score = score
        self.thisptr.setClassicScore(cpp_score)

    @property
    def borderAreaScore(PartialDecomposition self):
        """!@brief gets the border area score

        @note -1 iff not calculated yet, \see GCGconshdlrDecompCalcBorderAreaScore
        @returns border area score.
        """
        cdef double result = self.thisptr.getBorderAreaScore()
        return result


    @borderAreaScore.setter
    def borderAreaScore(PartialDecomposition self, double score):
        """!@brief set the border area score.
        """
        cdef double cpp_score = score
        self.thisptr.setBorderAreaScore(cpp_score)

    @property
    def maxWhiteScore(PartialDecomposition self):
        """!@brief gets the maximum white area score

        "maximum white score" is fraction of the area of the decomposed matrix that is neither block or border
        @note -1 iff not calculated yet, \see GCGconshdlrDecompCalcMaxWhiteScore
        @returns maximum  white area score
        .
        """
        cdef double result = self.thisptr.getMaxWhiteScore()
        return result


    @maxWhiteScore.setter
    def maxWhiteScore(PartialDecomposition self, double score):
        """!@brief set the maximum white area score.
        """
        cdef double cpp_score = score
        self.thisptr.setMaxWhiteScore(cpp_score)

    @property
    def maxForWhiteScore(PartialDecomposition self):
        """!@brief gets the maximum foreseeing white area score

        @note -1 iff not calculated yet, \see GCGconshdlrDecompCalcMaxForseeingWhiteScore
        @returns maximum foreseeing white area score
        .
        """
        cdef double result = self.thisptr.getMaxForWhiteScore()
        return result


    @maxForWhiteScore.setter
    def maxForWhiteScore(PartialDecomposition self, double score):
        """!@brief set the maximum foreseeing white area score.
        """
        cdef double cpp_score = score
        self.thisptr.setMaxForWhiteScore(cpp_score)

    @property
    def partForWhiteScore(PartialDecomposition self):
        """!@brief gets the setpartitioning maximum foreseeing white area score

        @note -1 iff not calculated yet, \see GGCGconshdlrDecompCalcSetPartForseeingWhiteScore
        @returns setpartitioning maximum foreseeing white area score
        .
        """
        cdef double result = self.thisptr.getSetPartForWhiteScore()
        return result


    @partForWhiteScore.setter
    def partForWhiteScore(PartialDecomposition self, double score):
        """!@brief set the setpartitioning maximum foreseeing white area score.
        """
        cdef double cpp_score = score
        self.thisptr.setSetPartForWhiteScore(cpp_score)

    @property
    def maxForWhiteAggScore(PartialDecomposition self):
        """!@brief gets the maximum foreseeing white area score with respect to aggregatable blocks

        @note -1 iff not calculated yet, \see GCGconshdlrDecompCalcMaxForeseeingWhiteAggScore
        @returns maximum foreseeing white area score with respect to aggregatable blocks
        .
        """
        cdef double result = self.thisptr.getMaxForWhiteAggScore()
        return result


    @maxForWhiteAggScore.setter
    def maxForWhiteAggScore(PartialDecomposition self, double score):
        """!@brief set the maximum foreseeing white area score with respect to aggregatable blocks.
        """
        cdef double cpp_score = score
        self.thisptr.setMaxForWhiteAggScore(cpp_score)

    @property
    def partForWhiteAggScore(PartialDecomposition self):
        """!@brief gets the setpartitioning maximum foreseeing white area score with respect to aggregateable

        @note -1 iff not calculated yet, \see GCGconshdlrDecompCalcSetPartForWhiteAggScore
        @returns setpartitioning maximum foreseeing white area score with respect to aggregateable.
        """
        cdef double result = self.thisptr.getSetPartForWhiteAggScore()
        return result


    @partForWhiteAggScore.setter
    def partForWhiteAggScore(PartialDecomposition self, double score):
        """!@brief set the setpartitioning maximum foreseeing white area score with respect to aggregateable.
        """
        cdef double cpp_score = score
        self.thisptr.setSetPartForWhiteAggScore(cpp_score)

    @property
    def bendersScore(PartialDecomposition self):
        """!@brief gets the benders score

        @note -1 iff not calculated yet, \see GCGconshdlrDecompCalcBendersScore
        @returns benders score.
        """
        cdef double result = self.thisptr.getBendersScore()
        return result


    @bendersScore.setter
    def bendersScore(PartialDecomposition self, double score):
        """!@brief set the benders score.
        """
        cdef double cpp_score = score
        self.thisptr.setBendersScore(cpp_score)

    @property
    def strongDecompScore(PartialDecomposition self):
        """!@brief gets the strong decomposition score

        @note -1 iff not calculated yet, \see GCGconshdlrDecompCalcStrongDecompositionScore
        @returns strong decomposition score.
        """
        cdef double result = self.thisptr.getStrongDecompScore()
        return result


    @strongDecompScore.setter
    def strongDecompScore(PartialDecomposition self, double score):
        """!@brief set the strong decomposition score.
        """
        cdef double cpp_score = score
        self.thisptr.setStrongDecompScore(cpp_score)

    def prepare(PartialDecomposition self):
        """!sorts the partialdec and calculates a its implicit assignments, hashvalue and evaluation

        @returns SCIP_OKAY if the result is consistent, SCIP_ERROR if there was an inconsistency.
        """
        self.thisptr.prepare()

    def aggInfoCalculated(PartialDecomposition self):
        """!@brief Checks if the aggregation information was already calculated
        @return true iff the aggregation information was already calculated.
        """
        cdef bool result = self.thisptr.aggInfoCalculated()
        return result


    def calcAggregationInformation(PartialDecomposition self, bool ignoreDetectionLimits):
        """!@brief computes if aggregation of sub problems is possible

        checks if aggregation of sub problems is possible and stores the corresponding aggregation information

        @param ignoreDetectionLimits Set to true if computation should ignore detection limits.

        This parameter is ignored if the patched bliss version is not present.
        """
        cdef bool cpp_ignoreDetectionLimits = ignoreDetectionLimits
        self.thisptr.calcAggregationInformation(cpp_ignoreDetectionLimits)

    def getConssForBlocks(PartialDecomposition self):
        cdef vector[vector[int] ] result = self.thisptr.getConssForBlocks()
        return result


    def getTranslatedpartialdecid(PartialDecomposition self):
        cdef int result = self.thisptr.getTranslatedpartialdecid()
        return result


    def setTranslatedpartialdecid(PartialDecomposition self, int decid):
        cdef int cpp_decid = decid
        self.thisptr.setTranslatedpartialdecid(cpp_decid)

    def buildDecChainString(PartialDecomposition self, buffer):
        """!@brief creates a detector chain short string for this partialdec, is built from detector chain.
        """
        c_buffer = str_conversion(buffer)
        self.thisptr.buildDecChainString(c_buffer)

    # END AUTOGENERATED BLOCK

    def matrix(PartialDecomposition self):
        cdef map[pair[int, int], double] mat = self.thisptr.writeNonzeros()
        return mat

    def matrixarray(PartialDecomposition self, obj=False, b=False):
        dictionary = self.matrix()
        X = list()
        Y = list()
        coeffs = list()
        for key, value in dictionary.items():
            if obj == True and b == True:
                X.append(key[0])
                Y.append(key[1])
                coeffs.append(value)
            elif b == True and obj == False:
                if key[0] != 0:
                    X.append(key[0])
                    Y.append(key[1])
                    coeffs.append(value)
            elif obj == True and b == False:
                if key[1] != self.getNVars():
                    X.append(key[0])
                    Y.append(key[1])
                    coeffs.append(value)
            else:
                if key[0] != 0 and key[1] != self.getNVars():
                    X.append(key[0])
                    Y.append(key[1])
                    coeffs.append(value)
        return X, Y, coeffs

    def visualize(PartialDecomposition self, fname=None, figsize=(12, 8), dpi=None, nonzero=True, only_boxes=False, obj=False, b=False, boxes=True, s=1, alpha=1, linkingcolor='#FFB72D', mastercolor='#1340C7', blockcolor='#718CDB', stairlinkingcolor='#886100', opencolor='#FFD88F', linecolor='#000000'):
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches

            plt.ioff()

            X, Y, coeffs = self.matrixarrayRhsAndObj(obj=obj, b=b)
            fig, ax = plt.subplots(figsize=figsize)

            #create the boxes
            if boxes==True:
                rowboxcounter = 0
                colboxcounter = 0

                if self.getNLinkingvars()!=0:
                    lvars = patches.Rectangle((0,0), self.getNLinkingvars(), self.getNConss(), linewidth=0.8, alpha=alpha, facecolor=linkingcolor, zorder=0, edgecolor=linecolor)
                    ax.add_patch(lvars)
                    colboxcounter+=self.getNLinkingvars()

                if self.getNMasterconss()!=0:
                    master = patches.Rectangle((0,0), self.getNVars(), self.getNMasterconss(), linewidth=0.8, alpha=alpha, facecolor=mastercolor, zorder=0, edgecolor=linecolor)
                    ax.add_patch(master)
                    rowboxcounter+=self.getNMasterconss()

                if self.getNMastervars()!=0:
                    colboxcounter+=self.getNMastervars()

                for b in range(self.getNBlocks()):
                    block = patches.Rectangle((colboxcounter,rowboxcounter), self.getNVarsForBlock(b), self.getNConssForBlock(b), linewidth=0.8, alpha=alpha, facecolor=blockcolor, zorder=0, edgecolor=linecolor)
                    ax.add_patch(block)
                    colboxcounter += self.getNVarsForBlock(b)
                    if self.getNStairlinkingvars(b)!=0:
                        stairlinking = patches.Rectangle((colboxcounter,rowboxcounter), self.getNStairlinkingvars(b), self.getNConssForBlock(b)+self.getNConssForBlock(b+1), linewidth=0.8, alpha=alpha, facecolor=stairlinkingcolor, zorder=0, edgecolor=linecolor)
                        ax.add_patch(stairlinking)
                    colboxcounter += self.getNStairlinkingvars(b)
                    rowboxcounter += self.getNConssForBlock(b)

                if self.getNOpenvars()!=0:
                    openrec = patches.Rectangle((colboxcounter,rowboxcounter), self.getNOpenvars(), self.getNOpenconss(), linewidth=0.8, alpha=alpha, facecolor=opencolor, zorder=0, edgecolor=linecolor)
                    ax.add_patch(openrec)
                    colboxcounter += self.getNOpenvars()
                    rowboxcounter += self.getNOpenconss()

            #set the zorder of the scatter-points
            zorderForBoxes = 0
            if boxes == True:
                zorderForBoxes = 1

            #plot the coefficients
            if only_boxes != True:
                if nonzero==True:
                    scatter=ax.scatter([y+0.5 for y in Y], [x+0.5 for x in X], c='black', s=s, alpha=1, zorder=zorderForBoxes)
                else:
                    if nonzero == False:
                        scatter=ax.scatter([y+0.5 for y in Y], [x+0.5 for x in X], c=coeffs, s=s, alpha=1, zorder=zorderForBoxes)
                    else:
                        scatter=ax.scatter([y+0.5 for y in Y], [x+0.5 for x in X], c=coeffs, cmap=nonzero, s=s, alpha=1, zorder=zorderForBoxes)
                    fig.colorbar(scatter)

            #adjust x-axis and y-axis
            if b == True:
                ax.set_xlim([0, self.getNVars()+1])
            else:
                ax.set_xlim([0, self.getNVars()])
            
            if obj == True:
                ax.set_ylim([-1, self.getNConss()])
            else:
                ax.set_ylim([0, self.getNConss()])
            ax.xaxis.tick_top()
            ax.invert_yaxis()

            if fname != None:
                plt.savefig(fname=fname, dpi=dpi)
                plt.close()
            else:
                plt.show()

        except ImportError:
            print("matplotlib is needed")

    def _repr_svg_(self):
        return self.__generate_visualization("svg")

    def _repr_png_(self):
        return self.__generate_visualization("png")

    cdef __generate_visualization(self, format="svg"):
        format = format.lower()
        if format not in ["svg", "png"]:
            raise ValueError(f"Format {format} is not supported. Only \"svg\" and \"png\" are supported.")

        if format not in self._visualizations:
            with tempfile.TemporaryDirectory() as td:
                temp_path = Path(td)

                gp_filename = temp_path.joinpath("vis.gp")
                outfile = temp_path.joinpath("vis").with_suffix(f".{format}")

                c_gp_filename = str_conversion(str(gp_filename))
                c_outfile = str_conversion(str(outfile))

                if format == "svg":
                    c_output_format = GP_OUTPUT_FORMAT_SVG
                elif format == "png":
                    c_output_format = GP_OUTPUT_FORMAT_PNG

                self.thisptr.generateVisualisation(c_gp_filename, c_outfile, c_output_format)

                if format == "svg":
                    data = outfile.read_text()
                elif format == "png":
                    data = outfile.read_bytes()
                self._visualizations[format] = data

        return self._visualizations[format]


    def __repr__(PartialDecomposition self):
        return f"<PartialDecomposition: nBlocks={self.getNBlocks()}, nMasterConss={self.getNMasterconss()}, nMasterVars={self.getNMastervars()}, nLinkingVars={self.getNLinkingvars()}, maxForWhiteScore={self.maxForWhiteScore}>"


cdef class DetProbData:
    """!class to manage the detection process and data for one coefficient matrix of a MIP, usually there is one detprobdata for the original and one detprobdata for the presolved problem.
    """
    cdef DETPROBDATA * thisptr
    cdef bool delete_thisptr

    # make DetProbData weak referentiable
    cdef object __weakref__

    def __cinit__(self):
        self.thisptr = NULL
        self.delete_thisptr = True

    def __dealloc__(self):
        if self.delete_thisptr and self.thisptr != NULL:
            del self.thisptr

    @staticmethod
    cdef create(DETPROBDATA* thisptr):
        if thisptr == NULL:
            raise Warning("cannot create DetProbData with DETPROBDATA* == NULL")
        new_DetProbData = DetProbData()
        new_DetProbData.thisptr = thisptr
        new_DetProbData.delete_thisptr = False
        return new_DetProbData

    # def __init__(DetProbData self, SCIP * scip, unsigned int _originalProblem):
    #     """@brief constructor
    #     @param scip SCIP data structure
    #     @param _originalProblem true iff the detprobdata is created for the presolved problem.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    @property
    def candidatesNBlocks(DetProbData self):
        """!< candidate for the number of blocks, second int indicates how often a candidate was added.
        """
        cdef vector[pair[int, int] ] result = self.thisptr.candidatesNBlocks
        return result

    
    @candidatesNBlocks.setter
    def candidatesNBlocks(DetProbData self, object candidatesNBlocks):
        cdef vector[pair[int, int] ] cpp_candidatesNBlocks = candidatesNBlocks
        self.thisptr.candidatesNBlocks = cpp_candidatesNBlocks

    @property
    def conspartitioncollection(DetProbData self):
        """!< collection of different constraint class distributions.
        """
        cdef vector[ConsPartition *] result = self.thisptr.conspartitioncollection
        return [ConsPart.create(r, <DetProbData>weakref.proxy(self)) for r in result]
    
    @conspartitioncollection.setter
    def conspartitioncollection(DetProbData self, object conspartitioncollection):
        # this seems to be possible only when we use C++11 (-std=c++11)
        # maybe it will be fixed in a future version of Cython
        cdef vector[ConsPartition *] cpp_conspartitioncollection
        cdef ConsPartition * conspartitioncollection_ptr = NULL
        cdef ConsPart conspartitioncollection_element
        for conspartitioncollection_element in conspartitioncollection:
            conspartitioncollection_ptr = <ConsPartition*> conspartitioncollection_element.thisptr
            cpp_conspartitioncollection.push_back(conspartitioncollection_ptr)
        self.thisptr.conspartitioncollection = cpp_conspartitioncollection

    @property
    def varpartitioncollection(DetProbData self):
        """!< collection of different variable class distributions.
        """
        cdef vector[VarPartition *] result = self.thisptr.varpartitioncollection
        return [VarPart.create(r) for r in result]
    
    @varpartitioncollection.setter
    def varpartitioncollection(DetProbData self, object varpartitioncollection):
        # this seems to be possible only when we use C++11 (-std=c++11)
        # maybe it will be fixed in a future version of Cython
        cdef vector[VarPartition *] cpp_varpartitioncollection
        cdef VarPartition * varpartitioncollection_ptr = NULL
        cdef VarPart varpartitioncollection_element
        for varpartitioncollection_element in varpartitioncollection:
            varpartitioncollection_ptr = <VarPartition*> varpartitioncollection_element.thisptr
            cpp_varpartitioncollection.push_back(varpartitioncollection_ptr)
        self.thisptr.varpartitioncollection = cpp_varpartitioncollection

    @property
    def classificationtime(DetProbData self):
        """!< time that was consumed by the classification of the constraint and variables classifiers.
        """
        cdef double result = self.thisptr.classificationtime
        return result

    
    @classificationtime.setter
    def classificationtime(DetProbData self, double classificationtime):
        cdef double cpp_classificationtime = classificationtime
        self.thisptr.classificationtime = cpp_classificationtime

    @property
    def nblockscandidatescalctime(DetProbData self):
        """!< time that was used to calulate the candidates of te block number.
        """
        cdef double result = self.thisptr.nblockscandidatescalctime
        return result

    
    @nblockscandidatescalctime.setter
    def nblockscandidatescalctime(DetProbData self, double nblockscandidatescalctime):
        cdef double cpp_nblockscandidatescalctime = nblockscandidatescalctime
        self.thisptr.nblockscandidatescalctime = cpp_nblockscandidatescalctime

    @property
    def postprocessingtime(DetProbData self):
        """!< time that was spent in postproceesing decomposigtions.
        """
        cdef double result = self.thisptr.postprocessingtime
        return result

    
    @postprocessingtime.setter
    def postprocessingtime(DetProbData self, double postprocessingtime):
        cdef double cpp_postprocessingtime = postprocessingtime
        self.thisptr.postprocessingtime = cpp_postprocessingtime

    @property
    def translatingtime(DetProbData self):
        """!< time that was spent by transforming partialdecs between presolved and orig problem.
        """
        cdef double result = self.thisptr.translatingtime
        return result

    
    @translatingtime.setter
    def translatingtime(DetProbData self, double translatingtime):
        cdef double cpp_translatingtime = translatingtime
        self.thisptr.translatingtime = cpp_translatingtime


    def addConsPartition(DetProbData self, ConsPart partition):
        """!@brief adds a constraint partition if it is no duplicate of an existing constraint partition.
        """
        cdef ConsPartition * cpp_partition = partition.thisptr
        self.thisptr.addConsPartition(cpp_partition)

    def addCandidatesNBlocksNVotes(DetProbData self, int candidate, int nvotes):
        """!@brief adds a candidate for block number and counts how often a candidate is added.
        """
        cdef int cpp_candidate = candidate
        cdef int cpp_nvotes = nvotes
        self.thisptr.addCandidatesNBlocksNVotes(cpp_candidate, cpp_nvotes)

    def addPartialdecToAncestor(DetProbData self, PartialDecomposition partialdec):
        """!@brief adds a partialdec to ancestor partialdecs
        @param partialdec partialdec that is added to the ancestor partialdecs.
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        self.thisptr.addPartialdecToAncestor(cpp_partialdec)

    def addPartialdecToOpen(DetProbData self, PartialDecomposition partialdec):
        """!@brief adds a partialdec to current partialdecs (data structure for partialdecs that are goin to processed in the propagation rounds)
        @param partialdec pointer of partialdec to be added
        @returns TRUE if the partialdecs was successfully added (i.

        e. it is no duplicate of a known partialdec)
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        cdef bool result = self.thisptr.addPartialdecToOpen(cpp_partialdec)
        return result


    def addPartialdecToFinished(DetProbData self, PartialDecomposition partialdec):
        """!@brief adds a partialdec to finished partialdecs
        @param partialdec pointer of partialdec that is going to be added to the finished partialdecs (data structure to carry finished decompositions)
        @returns TRUE if the partialdecs was successfully added (i.

        e. it is no duplicate of a known partialdec)
        @see addPartialdecToFinishedUnchecked()
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        cdef bool result = self.thisptr.addPartialdecToFinished(cpp_partialdec)
        return result


    def addPartialdecToFinishedUnchecked(DetProbData self, PartialDecomposition partialdec):
        """!@brief adds a partialdec to finished partialdecs without checking for duplicates, dev has to check this on his own
        @param partialdec pointer of partialdec that is going to be added unchecked to the finished partialdecs (data structure to carry finished decompositions)
        @see addPartialdecToFinished().
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        self.thisptr.addPartialdecToFinishedUnchecked(cpp_partialdec)

    def addVarPartition(DetProbData self, VarPart partition):
        """!@brief  adds a variable partition if it is no duplicate of an existing variable partition
        @param partition varpartition to be added.
        """
        cdef VarPartition * cpp_partition = partition.thisptr
        self.thisptr.addVarPartition(cpp_partition)

    def clearAncestorPartialdecs(DetProbData self):
        """!@brief clears ancestor partialdec data structure,
        @note does not free the partialdecs themselves.
        """
        self.thisptr.clearAncestorPartialdecs()

    def clearCurrentPartialdecs(DetProbData self):
        """!@brief clears current partialdec data structure
        @note does not free the partialdecs themselves.
        """
        self.thisptr.clearCurrentPartialdecs()

    def clearFinishedPartialdecs(DetProbData self):
        """!@brief clears finished partialdec data structure

        @note does not free the partialdecs themselves.
        """
        self.thisptr.clearFinishedPartialdecs()

    def createConssAdjacency(DetProbData self):
        """!@brief create the constraint adjacency datastructure that is used (if created) for some methods to faster access the constarints that have variables in common.
        """
        self.thisptr.createConssAdjacency()

    def freeTemporaryData(DetProbData self):
        """!@brief frees temporary data that is only needed during the detection process.
        """
        self.thisptr.freeTemporaryData()

    def getAncestorPartialdec(DetProbData self, int partialdecindex):
        """!@brief returns a partialdec from ancestor partialdec data structure with given index

        @returns partialdec from ancestor partialdec data structure.
        """
        cdef int cpp_partialdecindex = partialdecindex
        cdef PARTIALDECOMP * result = self.thisptr.getAncestorPartialdec(cpp_partialdecindex)
        return PartialDecomposition.create(result)

    def getConsPartition(DetProbData self, int partitionIndex):
        """!@brief returns pointer to a constraint partition
        @return pointer to a cosntraint partition with the given index.
        """
        cdef int cpp_partitionIndex = partitionIndex
        cdef ConsPartition * result = self.thisptr.getConsPartition(cpp_partitionIndex)
        return ConsPart.create(result, <DetProbData>weakref.proxy(self))

    def getCons(DetProbData self, int consIndex):
        """!@brief returns the SCIP constraint related to a constraint index
        @return the SCIP constraint related to a constraint index.
        """
        return Constraint.create(self.thisptr.getCons(consIndex))

    def getConssForCons(DetProbData self, int consIndex):
        """!@brief return array of constraint indices that have a common variable with the given constraint
        @return return vector of constraint indices that have a common variable with the given constraint
        @note constraint adjacency data structure has to initilized.
        """
        cdef int cpp_consIndex = consIndex
        cdef vector[int] result = self.thisptr.getConssForCons(cpp_consIndex)
        return result


    def getConssForVar(DetProbData self, int varIndex):
        """!\brief returns the constraint indices of the coefficient matrix for a variable
        @return vector of constraint indices that have a nonzero entry with this variable.
        """
        cdef int cpp_varIndex = varIndex
        cdef vector[int] result = self.thisptr.getConssForVar(cpp_varIndex)
        return result


    def getOpenPartialdecs(DetProbData self):
        """!@brief determines all partialdecs from current (open) partialdec data structure
        @returns  all partialdecs in current (open) partialdec data structure
        /
        /**.
        """
        cdef vector[PARTIALDECOMP *] result = self.thisptr.getOpenPartialdecs()
        return [PartialDecomposition.create(r) for r in result]

    def getFinishedPartialdec(DetProbData self, int partialdecindex):
        """!@brief returns a partialdec from finished partialdec data structure
        @return  partialdec from finished partialdec data structure.
        """
        cdef int cpp_partialdecindex = partialdecindex
        cdef PARTIALDECOMP * result = self.thisptr.getFinishedPartialdec(cpp_partialdecindex)
        return PartialDecomposition.create(result)

    def getFinishedPartialdecs(DetProbData self):
        """!@brief gets all finished partialdecs
        @returns all finished partialdecs.
        """
        cdef vector[PARTIALDECOMP *] result = self.thisptr.getFinishedPartialdecs()
        return [PartialDecomposition.create(r) for r in result]

    def getIndexForCons(DetProbData self, Constraint cons):
        """@brief returns the constraint index related to a SCIP constraint
        @param cons the SCIP constraint pointer the index is asked for
        @return the constraint index related to a SCIP constraint.
        """
        return self.thisptr.getIndexForCons(cons.scip_cons)

    # def getIndexForVar(DetProbData self, SCIP_VAR * var):
    #     """@brief returns the variable index related to a SCIP variable
    #     @param var variable pointer the index is asked for
    #     @return the variable index related to a SCIP variable.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def getNAncestorPartialdecs(DetProbData self):
        """!@brief returns size of ancestor partialdec data structure
        @return size of ancestor partialdec data structure.
        """
        cdef int result = self.thisptr.getNAncestorPartialdecs()
        return result


    def getNConsPartitions(DetProbData self):
        """!@brief returns number of different constraint partitions
        @return number of different constraint partitions.
        """
        cdef int result = self.thisptr.getNConsPartitions()
        return result


    def getNConss(DetProbData self):
        """!@brief returns the number of variables considered in the detprobdata
        @return number of variables considered in the detprobdata.
        """
        cdef int result = self.thisptr.getNConss()
        return result


    def getNConssForCons(DetProbData self, int consIndex):
        """!@brief returns the number of constraints for a given constraint
        @return the number of constraints for a given constraint.
        """
        cdef int cpp_consIndex = consIndex
        cdef int result = self.thisptr.getNConssForCons(cpp_consIndex)
        return result


    def getNConssForVar(DetProbData self, int varIndex):
        """!@brief returns the number of constraints for a given variable where the var has a nonzero entry in
        @return the number of constraints for a given variable.
        """
        cdef int cpp_varIndex = varIndex
        cdef int result = self.thisptr.getNConssForVar(cpp_varIndex)
        return result


    def getNOpenPartialdecs(DetProbData self):
        """!@brief returns size of current (open) partialdec data structure
        @return size of current (open) partialdec data structure.
        """
        cdef int result = self.thisptr.getNOpenPartialdecs()
        return result


    def getNFinishedPartialdecs(DetProbData self):
        """!returns size of finished partialdec data structure
        @return  size of finished partialdec data structure.
        """
        cdef int result = self.thisptr.getNFinishedPartialdecs()
        return result


    def getNPartialdecs(DetProbData self):
        """!returns the number of stored partialdecs
        @return  number of stored partialdecs.
        """
        cdef int result = self.thisptr.getNPartialdecs()
        return result


    def getNNonzeros(DetProbData self):
        """!@brief returns the number of nonzero entries in the coefficient matrix
        @return the number of nonzero entries in the coefficient matrix.
        """
        cdef int result = self.thisptr.getNNonzeros()
        return result


    def getNVarPartitions(DetProbData self):
        """!@brief returns number of different variable partitions
        @return  number of different variable partitions.
        """
        cdef int result = self.thisptr.getNVarPartitions()
        return result


    def getNVars(DetProbData self):
        """!@brief return the number of variables considered in the detprobdata
        @return the number of variables considered in the detprobdata.
        """
        cdef int result = self.thisptr.getNVars()
        return result


    def getNVarsForCons(DetProbData self, int consIndex):
        """!@brief returns the number of variables for a given constraint
        @return the number of variables for a given constraint.
        """
        cdef int cpp_consIndex = consIndex
        cdef int result = self.thisptr.getNVarsForCons(cpp_consIndex)
        return result


    def getOrigVarsFixedZero(DetProbData self):
        """!@brief returns pointers to all orig vars that are fixed to zero
        @returns vector of vars.
        """
        cdef vector[SCIP_VAR *] result = self.thisptr.getOrigVarsFixedZero()
        return [Variable.create(v) for v in result]


    def getRelevantConss(DetProbData self):
        """!@brief returns pointers to all constraints that are not marked as deleted or obsolete
        @returns vector of conss.
        """
        cdef vector[SCIP_CONS *] result = self.thisptr.getRelevantConss()
        return [Constraint.create(c) for c in result]


    def getRelevantVars(DetProbData self):
        """!@brief returns pointers to all problem vars that are not fixed to 0
        @returns vector of vars.
        """
        cdef vector[SCIP_VAR *] result = self.thisptr.getRelevantVars()
        return [Variable.create(v) for v in result]



    def getModel(DetProbData self):
        """!@brief returns the corresponding Model instance wrapping the scip data structure
        @return the corresponding Model instance wrapping scip data structure.
        """
        cdef SCIP * scip = self.thisptr.getScip()
        return Model.create(scip)

    def getSortedCandidatesNBlocks(DetProbData self, object candidates):
        """!@brief gets the candidates for number of blocks added by the user followed by the found ones sorted in descending order by how often a candidate was proposed
        @param candidates will contain the candidates for number of blocks sorted in descending order by how often a candidate was added.
        """
        cdef vector[int] cpp_candidates = candidates
        self.thisptr.getSortedCandidatesNBlocks(cpp_candidates)

    def getVal(DetProbData self, int row, int col):
        """!@brief returns a coefficient from the coefficient matrix
        @return a coefficient from the coefficient matrix.
        """
        cdef int cpp_row = row
        cdef int cpp_col = col
        cdef double result = self.thisptr.getVal(cpp_row, cpp_col)
        return result


    def getValsForCons(DetProbData self, int consIndex):
        """!@brief returns the nonzero coefficients of the coefficient matrix for a constraint
        @return vector of coefficients of in matrix for constraints
        @note same order as in @see getVarsForCons().
        """
        cdef int cpp_consIndex = consIndex
        cdef vector[double] result = self.thisptr.getValsForCons(cpp_consIndex)
        return result


    def getVarPartition(DetProbData self, int partitionIndex):
        """!@brief returns pointer to a variable partition with given index
        @return pointer to a variable partition with given index.
        """
        cdef int cpp_partitionIndex = partitionIndex
        cdef VarPartition * result = self.thisptr.getVarPartition(cpp_partitionIndex)
        return VarPart.create(result)

    def getVarPartitions(DetProbData self):
        """!@brief returns vector to stored variable partitions
        @return returns vector to stored variable partitions.
        """
        cdef vector[VarPartition *] result = self.thisptr.getVarPartitions()
        return [VarPart.create(r) for r in result]

    def getVar(DetProbData self, int varIndex):
        """!@brief returns SCIP variable related to a variable index
        @return SCIP variable pointer related to a variable index.
        """
        cdef int cpp_varIndex = varIndex
        # TODO implement function
        raise NotImplementedError()

    def getVarsForCons(DetProbData self, int consIndex):
        """!@brief returns the variable indices of the coefficient matrix for a constraint
        @return the variable indices of the coefficient matrix for a constraint.
        """
        cdef int cpp_consIndex = consIndex
        cdef vector[int] result = self.thisptr.getVarsForCons(cpp_consIndex)
        return result


    def isConsCardinalityCons(DetProbData self, int consindexd):
        """!@brief returns whether a constraint is a cardinality constraint, i.

        e. of the \f$\sum_{i} x_i = b\f$
        @param consindexd index of constraint that is be checked
        @return returns whether a constraint is a cardinality constraint
        """
        cdef int cpp_consindexd = consindexd
        cdef bool result = self.thisptr.isConsCardinalityCons(cpp_consindexd)
        return result


    def isConssAdjInitialized(DetProbData self):
        """!@brief determines whether or not the constraint-constraint adjacency data structure is initilized

        @returns true iff the constraint-constraint adjacency data structure is initilized.
        """
        cdef unsigned int result = self.thisptr.isConssAdjInitialized()
        return result


    def isConsSetpp(DetProbData self, int consindexd):
        """!@brief is cons with specified indec partitioning, or packing covering constraint?
        @param consindexd index of the given cons
        @return is cons with specified indec partitioning, or packing covering constraint.
        """
        cdef int cpp_consindexd = consindexd
        cdef bool result = self.thisptr.isConsSetpp(cpp_consindexd)
        return result


    def isConsSetppc(DetProbData self, int consindexd):
        """!@brief is cons with specified index partitioning packing, or covering constraint?
        @param consindexd index of cons to be checked
        @return whether a constraint is partitioning packing, or covering constraint?.
        """
        cdef int cpp_consindexd = consindexd
        cdef bool result = self.thisptr.isConsSetppc(cpp_consindexd)
        return result


    # def isFiniteNonnegativeIntegral(DetProbData self, SCIP * scip, double x):
    #     """@brief is constraint ranged row, i.

    #     e., -inf < lhs < rhs < inf?

    #     @returns whether val is ranged row
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def isPartialdecDuplicateofFinished(DetProbData self, PartialDecomposition partialdec):
        """!@brief check if partialdec is a duplicate of an existing finished partialdec
        @param partialdec partialdec to be checked
        @returns TRUE iff partialdec is a duplicate of an existing finished partialdec.
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        cdef unsigned int result = self.thisptr.isPartialdecDuplicateofFinished(cpp_partialdec)
        return result


    def isAssignedToOrigProb(DetProbData self):
        """!@brief returns true if the matrix structure corresponds to the presolved problem
        @return TRUE if the matrix structure corresponds to the presolved problem.
        """
        cdef unsigned int result = self.thisptr.isAssignedToOrigProb()
        return result


    # def isRangedRow(DetProbData self, SCIP * scip, double lhs, double rhs):
    #     """@brief is constraint ranged row, i.

    #     e., -inf < lhs < rhs < inf?
    #     @returns whether constraint is ranged row
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def partialdecIsNoDuplicateOfPartialdecs(DetProbData self, PartialDecomposition comppartialdec, object partialdecs, bool sort):
        """!@brief check if partialdec is a duplicate of any given partialdecs
        @param comppartialdec partialdec to be checked
        @param partialdecs partialdecs to compare with
        @param sort sort the vars and conss data structures in the partialdecs by their indices
        @return TRUE iff partialdec is no duplicate of any given partialdecs.
        """
        cdef PARTIALDECOMP * cpp_comppartialdec = comppartialdec.thisptr
        # this seems to be possible only when we use C++11 (-std=c++11)
        # maybe it will be fixed in a future version of Cython
        cdef vector[PARTIALDECOMP *] cpp_partialdecs
        cdef PARTIALDECOMP * partialdecs_ptr = NULL
        cdef PartialDecomposition partialdecs_element
        for partialdecs_element in partialdecs:
            partialdecs_ptr = <PARTIALDECOMP*> partialdecs_element.thisptr
            cpp_partialdecs.push_back(partialdecs_ptr)
        cdef bool cpp_sort = sort
        cdef unsigned int result = self.thisptr.partialdecIsNoDuplicateOfPartialdecs(cpp_comppartialdec, cpp_partialdecs, cpp_sort)
        return result


    # def printBlockcandidateInformation(DetProbData self, SCIP * scip, FILE * file):
    #     """@brief output method for json file writer to write block candidate information
    #     @param scip SCIP data structure
    #     @param file  output file or NULL for standard output.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def printPartitionInformation(DetProbData self, FILE * file):
    #     """@brief output method for json file writer to write partition candidate information
    #     @param file output file or NULL for standard output.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def sortFinishedForScore(DetProbData self):
        """!@brief sorts partialdecs in finished partialdecs data structure according to the current scoretype.
        """
        self.thisptr.sortFinishedForScore()

    def translatePartialdecs(DetProbData self, DetProbData otherdata, object otherpartialdecs):
        """!@brief translates partialdecs if the index structure of the problem has changed, e.

        g. due to presolving
        @return translated partialdecs
        """
        cdef DETPROBDATA * cpp_otherdata = otherdata.thisptr
        # this seems to be possible only when we use C++11 (-std=c++11)
        # maybe it will be fixed in a future version of Cython
        cdef vector[PARTIALDECOMP *] cpp_otherpartialdecs
        cdef PARTIALDECOMP * otherpartialdecs_ptr = NULL
        cdef PartialDecomposition otherpartialdecs_element
        for otherpartialdecs_element in otherpartialdecs:
            otherpartialdecs_ptr = <PARTIALDECOMP*> otherpartialdecs_element.thisptr
            cpp_otherpartialdecs.push_back(otherpartialdecs_ptr)
        cdef vector[PARTIALDECOMP *] result = self.thisptr.translatePartialdecs(cpp_otherdata, cpp_otherpartialdecs)
        return [PartialDecomposition.create(r) for r in result]


cdef class ConsPart:
    cdef ConsPartition * thisptr
    cdef DetProbData detProbData
    cdef bool delete_thisptr

    def __cinit__(self):
        self.thisptr = NULL
        self.delete_thisptr = True

    def __dealloc__(self):
        if self.delete_thisptr and self.thisptr != NULL:
            del self.thisptr

    @staticmethod
    cdef create(ConsPartition* thisptr, DetProbData detProbData):
        if thisptr == NULL:
            raise Warning("cannot create ConsPart with ConsPartition* == NULL")
        new_ConsPart = ConsPart()
        new_ConsPart.thisptr = thisptr
        new_ConsPart.delete_thisptr = False
        new_ConsPart.detProbData = detProbData
        return new_ConsPart


    # def addClass(ConsPart self, name, desc, CONS_DECOMPINFO decompInfo):
    #     """creates a new class, returns index of the class.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def assignConsToClass(ConsPart self, Constraint cons, int classindex):
        """!assigns a constraint to a class.
        """
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        cdef int cpp_classindex = classindex
        self.thisptr.assignConsToClass(cpp_consindex, cpp_classindex)

    def getAllSubsets(ConsPart self, bool both, bool only_master, bool only_pricing):
        """!returns a vector containing all possible subsets of the chosen classindices.
        """
        cdef bool cpp_both = both
        cdef bool cpp_only_master = only_master
        cdef bool cpp_only_pricing = only_pricing
        cdef vector[vector[int] ] result = self.thisptr.getAllSubsets(cpp_both, cpp_only_master, cpp_only_pricing)
        return result


    # def getClassDecompInfo(ConsPart self, int classindex):
    #     """returns the decomposition info of a class.
    #     """
    #     cdef int cpp_classindex = classindex
    #     # TODO implement function
    #     raise NotImplementedError()

    def getClassNameOfCons(ConsPart self, Constraint cons):
        """!returns the name of the class a constraint is assigned to.
        """
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        cdef const char * result = self.thisptr.getClassNameOfCons(cpp_consindex)
        return result.decode('utf-8')


    def getClassOfCons(ConsPart self, Constraint cons):
        """!returns the index of the class a constraint is assigned to.
        """
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        cdef int result = self.thisptr.getClassOfCons(cpp_consindex)
        return result


    # def getConssToClasses(ConsPart self):
    #     """returns vector containing the assigned class of each constraint.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def getNConss(ConsPart self):
        """!returns the number of constraints.
        """
        cdef int result = self.thisptr.getNConss()
        return result


    def getNConssOfClasses(ConsPart self):
        """!returns a vector with the numbers of constraints that are assigned to the classes.
        """
        cdef vector[int] result = self.thisptr.getNConssOfClasses()
        return result


    def isConsClassified(ConsPart self, Constraint cons):
        """!returns whether a constraint is already assigned to a class.
        """
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        cdef bool result = self.thisptr.isConsClassified(cpp_consindex)
        return result


    def reduceClasses(ConsPart self, int maxNumberOfClasses):
        """!returns partition with reduced number of classes
        if the current number of classes is greater than an upper bound
        and lower than 2*(upper bound) (returns NULL otherwise).
        """
        cdef int cpp_maxNumberOfClasses = maxNumberOfClasses
        cdef ConsPartition * result = self.thisptr.reduceClasses(cpp_maxNumberOfClasses)
        return ConsPart.create(result, self.detProbData)


    def getName(ConsPart self):
        """!returns the name of the partition"""
        return self.thisptr.getName().decode('utf-8')

    # def setClassDecompInfo(ConsPart self, int classindex, CONS_DECOMPINFO decompInfo):
    #     """sets the decomposition code of a class.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def __repr__(ConsPart self):
        return f"<ConsPart: name={self.getName()}>"


cdef class VarPart:
    cdef VarPartition * thisptr
    cdef bool delete_thisptr

    def __cinit__(self):
        self.thisptr = NULL
        self.delete_thisptr = True

    def __dealloc__(self):
        if self.delete_thisptr and self.thisptr != NULL:
            del self.thisptr

    @staticmethod
    cdef create(VarPartition* thisptr):
        if thisptr == NULL:
            raise Warning("cannot create VarPart with VarPartition* == NULL")
        new_VarPart = VarPart()
        new_VarPart.thisptr = thisptr
        new_VarPart.delete_thisptr = False
        return new_VarPart


    # def addClass(VarPart self, name, desc, VAR_DECOMPINFO decompInfo):
    #     """creates a new class, returns index of the class.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def assignVarToClass(VarPart self, int varindex, int classindex):
        """!assigns a variable to a class.
        """
        cdef int cpp_varindex = varindex
        cdef int cpp_classindex = classindex
        self.thisptr.assignVarToClass(cpp_varindex, cpp_classindex)

    def getAllSubsets(VarPart self, bool all, bool linking, bool master, bool block):
        """!returns a vector containing all possible subsets of the chosen classindices.
        """
        cdef bool cpp_all = all
        cdef bool cpp_linking = linking
        cdef bool cpp_master = master
        cdef bool cpp_block = block
        cdef vector[vector[int] ] result = self.thisptr.getAllSubsets(cpp_all, cpp_linking, cpp_master, cpp_block)
        return result


    # def getClassDecompInfo(VarPart self, int classindex):
    #     """returns the decomposition info of a class.
    #     """
    #     cdef int cpp_classindex = classindex
    #     # TODO implement function
    #     raise NotImplementedError()

    def getClassNameOfVar(VarPart self, int varindex):
        """!returns the name of the class a variable is assigned to.
        """
        cdef int cpp_varindex = varindex
        cdef const char * result = self.thisptr.getClassNameOfVar(cpp_varindex)
        return result.decode('utf-8')


    def getClassOfVar(VarPart self, int varindex):
        """!returns the index of the class a variable is assigned to.
        """
        cdef int cpp_varindex = varindex
        cdef int result = self.thisptr.getClassOfVar(cpp_varindex)
        return result


    # def getVarsToClasses(VarPart self):
    #     """returns vector containing the assigned class of each variable.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def getNVars(VarPart self):
        """!returns the number of variables.
        """
        cdef int result = self.thisptr.getNVars()
        return result


    def getNVarsOfClasses(VarPart self):
        """!returns a vector with the numbers of variables that are assigned to the classes.
        """
        cdef vector[int] result = self.thisptr.getNVarsOfClasses()
        return result


    def isVarClassified(VarPart self, int varindex):
        """!returns whether a variable is already assigned to a class.
        """
        cdef int cpp_varindex = varindex
        cdef bool result = self.thisptr.isVarClassified(cpp_varindex)
        return result


    def reduceClasses(VarPart self, int maxNumberOfClasses):
        """!returns partition with reduced number of classes
        if the current number of classes is greater than an upper bound
        and lower than 2*(upper bound) (returns NULL otherwise).
        """
        cdef int cpp_maxNumberOfClasses = maxNumberOfClasses
        cdef VarPartition * result = self.thisptr.reduceClasses(cpp_maxNumberOfClasses)
        return VarPart.create(result)

    # def setClassDecompInfo(VarPart self, int classindex, VAR_DECOMPINFO decompInfo):
    #     """sets the decomposition code of a class.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()
