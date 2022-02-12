cdef class PartialDecomposition:
    """class to manage partial decompositions

    each partialdec corresponds to one :class:`DetProbData` which contains the problem information,
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
    #     """Standard constructor, creates empty partialdec with unique id

    #     :note: initially, all conss and vars are open.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def __init__(PartialDecomposition self, PartialDecomposition partialdecToCopy):
    #     """copy constructor.

    #     """
    #     cdef PARTIALDECOMP * cpp_partialdecToCopy = partialdecToCopy.thisptr
    #     self.thisptr = new PARTIALDECOMP(cpp_partialdecToCopy)

    def __copy__(self):
        return PartialDecomposition.create(new PARTIALDECOMP(self.thisptr))

    def copy(self):
        return copy(self)

    def fixConsToMaster(PartialDecomposition self, Constraint cons):
        """Fixes a Constraint to the master constraints.

        :param cons: scip#Constraint to add
        :type cons: Constraint
        """
        self.thisptr.fixConsToMaster(cons.scip_cons)

    def fixConssToMaster(self, conss):
        """Fixes all constraints to the master constraints.

        :param conss: An iterable of scip#Constraint objects
        :type conss: An iterable of scip#Constraint objects
        :raises TypeError: occurs if conss is not an Iterable
        .. seealso:: * :meth:`fixConsToMaster`
        """
        if not isinstance(conss, Iterable):
            raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(conss)))
        for cons in conss:
            self.fixConsToMaster(<Constraint?>cons)

    def fixConsToBlockId(PartialDecomposition self, Constraint cons, int block_id):
        """Adds a constraint to a block.

        :param cons: scip#Constraint to add
        :type cons: Constraint
        :param block_id: id of block to add.
        :type block_id: int
        :note: The passed ``block_id`` has to by of type integer and use the internal numbering of blocks. Before
        calling this method, one has to ensure that the specified block exists. A more convenient alternative is
        fixConsToBlock().

        .. seealso:: * :meth:`addBlock`
        """
        self.thisptr.fixConsToBlock(cons.scip_cons, block_id)

    def fixConssToBlockId(self, conss, int block_id):
        """Adds all constraints to a block.

        :param conss: An iterable of scip#Constraint objects
        :param block_id: id of block to add.
        .. seealso:: * :meth:`fixConsToBlockId`
                     * :meth:`fixConssToBlock`
        """
        if not isinstance(conss, Iterable):
            raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(conss)))
        for cons in conss:
            self.fixConsToBlockId(<Constraint?>cons, block_id)

    def fixConsToBlock(self, Constraint cons, object block):
        """Adds a constraint to a block.

        :param cons: scip#Constraint to add
        :param block: identifier of block to add. Can be any hashable Python object.
        :return: the internal block_id assigned to the block.

        This method will automatically manage the blocks of the decomposition and create blocks if neccessary.
        The passed ``block`` identifiers will be mapped to internal block ids.

        To address the internal blocks directly, use :meth:`fixConsToBlockId`.

        .. seealso:: * :meth:`fixConssToBlock`
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
        """Adds all constraints to a block.

        :param conss: An iterable of scip#Constraint objects
        :param block: identifier of block to add. Can be any hashable Python object.
        :return: the internal block_id assigned to the block.

        .. seealso:: * :meth:`fixConsToBlock`
        """
        if not isinstance(conss, Iterable):
            raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(conss)))
        cdef int block_id = -1  #assign no block_id at first
        for cons in conss:
            # block_id is the same for all blocks
            block_id = self.fixConsToBlock(<Constraint?>cons, block)
        return block_id

    def getOpenconss(PartialDecomposition self):
        """Gets a vector containing constraint ids not assigned yet as vector

        :return: returns a vector containing constraint ids not assigned yet as vector.
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
        """Adds a block

        :return: the number (id) of the new block
        """
        cdef int result = self.thisptr.addBlock()
        return result

    def addClockTime(PartialDecomposition self, double clocktime):
        """Adds detection time of one detector

        incorporates the needed time of some detector in the detector chain.
        """
        cdef double cpp_clocktime = clocktime
        self.thisptr.addClockTime(cpp_clocktime)

    def addDecChangesFromAncestor(PartialDecomposition self, PartialDecomposition ancestor):
        """Adds the statistical differences to an ancestor

        incorporates the changes from ancestor partialdec into the statistical data structures.
        """
        cdef PARTIALDECOMP * cpp_ancestor = ancestor.thisptr
        self.thisptr.addDecChangesFromAncestor(cpp_ancestor)

    def addDetectorChainInfo(PartialDecomposition self, decinfo):
        """Add information about the detector chain

        adds a detectorchain information string to the corresponding vector
        (that carries information for each detector call)
        """
        c_decinfo = str_conversion(decinfo)
        self.thisptr.addDetectorChainInfo(c_decinfo)

    def addNNewBlocks(PartialDecomposition self, int nnewblocks):
        """Adds how many new blocks were introduced

        bookkeeping information: adds number of new blocks created by a detector added to detector chain.
        """
        cdef int cpp_nnewblocks = nnewblocks
        self.thisptr.addNNewBlocks(cpp_nnewblocks)

    def addPctConssFromFree(PartialDecomposition self, double pct):
        """Adds percentage of closed constraints

        bookkeeping information: fraction of constraints that are not longer open for a detector added to detector chain.
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssFromFree(cpp_pct)

    def addPctConssToBlock(PartialDecomposition self, double pct):
        """Adds percentage of constraints assigned to blocks

        bookkeeping information: adds fraction of constraints assigned to a block for a detector added to detector chain
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssToBlock(cpp_pct)

    def addPctConssToBorder(PartialDecomposition self, double pct):
        """Adds percentage of constraints assigned to border

        bookkeeping information: adds fraction of constraints assigned to the border for a detector added to detector chain.
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssToBorder(cpp_pct)

    def addPctVarsFromFree(PartialDecomposition self, double pct):
        """Adds percentage of closed variables

        bookkeeping information: adds fraction of variables that are not longer open for a detector added to detector chain.
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsFromFree(cpp_pct)

    def addPctVarsToBlock(PartialDecomposition self, double pct):
        """Adds percentage of variables assigned to blocks

        bookkeeping information: adds fraction of variables assigned to a block for a detector added to detector chain
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsToBlock(cpp_pct)

    def addPctVarsToBorder(PartialDecomposition self, double pct):
        """Adds percentage of variables assigned to border

        bookkeeping information: adds fraction of variables assigned to the border for a detector added to detector chain.
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsToBorder(cpp_pct)

    def alreadyAssignedConssToBlocks(PartialDecomposition self):
        """Method to check if at least one constraint is assigned to some block

        :return: True iff at least one constraint is assigned to a block
        """
        cdef bool result = self.thisptr.alreadyAssignedConssToBlocks()
        return result

    # def assignBorderFromConstoblock(PartialDecomposition self, SCIP_HASHMAP * constoblock, int givenNBlocks):
    #     """assigns open conss to master

    #     assigns open constraints to master according to the cons assignment information given in constoblock hashmap
    #     :return: scip return code
    #     :note: for conss assigned to blocks according to constoblock there is no assignment \see assignPartialdecFromConstoblock
    #     :note: master assignment is indicated by assigning cons to index additionalNBlocks
    #     .
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def assignCurrentStairlinking(PartialDecomposition self):
        """Assigns open vars to stairlinking if appropriate

        assigns open vars to stairlinking if they can be found in exactly two consecutive blocks

        :return: True iff at least one stairlinkingvar was assigned.
        """
        cdef bool result = self.thisptr.assignCurrentStairlinking()
        return result

    def assignOpenConssToMaster(PartialDecomposition self):
        """Assigns open conss to master.
        """
        self.thisptr.assignOpenConssToMaster()

    # def assignPartialdecFromConstoblock(PartialDecomposition self, SCIP_HASHMAP * constoblock, int additionalNBlocks):
    #     """assigns conss structure according to given hashmap

    #     adds blocks and assigns open conss to a new block or to master
    #     according to the cons assignment information given in constoblock hashmap
    #     :return: scip return code
    #     \see assignPartialdecFromConstoblockVector()
    #     :note: master assignment is indicated by assigning cons to index additionalNBlocks
    #     .
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def assignPartialdecFromConstoblockVector(PartialDecomposition self, object constoblock, int additionalNBlocks):
        """assigns conss structure according to given vector

        adds blocks and assigns open conss to a new block or to master
        according to the cons assignment information given in constoblock vector

        :return: scip return code
        .. seealso::  * :meth:`assignPartialdecFromConstoblock()`
        .. note:: master is indicated by assigning cons to index additionalNBlocks.
        """
        cdef vector[int] cpp_constoblock = constoblock
        cdef int cpp_additionalNBlocks = additionalNBlocks
        # TODO implement function
        raise NotImplementedError()

    def assignSmallestComponentsButOneConssAdjacency(PartialDecomposition self):
        """computes components by connectedness of conss and vars

        computes components corresponding to connectedness of conss and vars
        and assigns them accordingly (all but one of largest components)

        strategy: assigns all conss same block if they are connected
        two constraints are adjacent if there is a common variable

        .. note:: this relies on the consadjacency structure of the detprobdata
        hence it cannot be applied in presence of linking variables.
        """
        self.thisptr.assignSmallestComponentsButOneConssAdjacency()

    def calcStairlinkingVars(PartialDecomposition self):
        """Reassigns linking vars to stairlinkingvars if possible

        potentially reorders blocks for making a maximum number of linking vars stairlinking
        if all vars that connect exactly two blocks have a staircase structure, all of them become stairlinkingvars
        otherwise, the stairlinking assignment is done greedily
        .. note:: precondition: partialdec does not have any stairlinking vars.
        """
        self.thisptr.calcStairlinkingVars()

    def checkAllConssAssigned(PartialDecomposition self):
        """Checks if all conss are assigned

        returns True iff all constraints are assigned and deletes the vector open conss if so

        :return: True iff all constraints are assigned
        """
        cdef bool result = self.thisptr.checkAllConssAssigned()
        return result

    def checkConsistency(PartialDecomposition self):
        """Checks whether the assignments in the partialdec are consistent

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

        :return: True iff the partialdec seems to be consistent
        """
        cdef bool result = self.thisptr.checkConsistency()
        return result

    def complete(PartialDecomposition self):
        """Assigns all open constraints and open variables trivially

        strategy: assigns all open conss and vars to blocks if they can be refined there, otherwise to the master

        :note: partialdecomps should usually be completed by a detector, only use this function if you know what you are doing.
        """
        self.thisptr.complete()

    def completeByConnected(PartialDecomposition self):
        """Assigns all open constraints and open variables

        strategy: assigns all conss and vars to the same block if they are connected,
        a cons and a var are adjacent if the var appears in the cons.
        """
        self.thisptr.completeByConnected()

    def completeByConnectedConssAdjacency(PartialDecomposition self):
        """Assigns all open constraints and open variables

        strategy: assigns all conss and vars to the same block if they are connected
        a cons and a var are adjacent if the var appears in the cons
        .. note:: this relies on the consadjacency structure of the detprobdata
        hence it cannot be applied in presence of linking variables.
        """
        self.thisptr.completeByConnectedConssAdjacency()

    def completeGreedily(PartialDecomposition self):
        """Assigns all open constraints and open variables

        strategy: assigns a cons (and related vars) to a new block if possible,
        if not to an existing block if possible (by means of prior var assignments)
        and finally to master, if there does not exist such a block.
        """
        self.thisptr.completeGreedily()

    def removeMastercons(PartialDecomposition self, int consid):
        """Removes the given cons from master.
        """
        cdef int cpp_consid = consid
        self.thisptr.removeMastercons(cpp_consid)

    def considerImplicits(PartialDecomposition self):
        """Assigns every open cons/var.

        Assignments happen as follows:
        - to the respective block if it hits exactly one blockvar/blockcons and no open vars/conss
        - to master/linking if it hits blockvars/blockconss assigned to different blocks
        - and every cons to master that hits a master var
        - and every var to master if it does not hit any blockcons and has no open cons
        - leave the cons/variableopen if nothing from the above holds
        """
        self.thisptr.considerImplicits()

    def copyPartitionStatistics(PartialDecomposition self, PartialDecomposition otherpartialdec):
        """Copies the given partialdec's partition statistics

        :param otherpartialdec: partialdec whose partition statistics are to be copied.
        """
        cdef PARTIALDECOMP * cpp_otherpartialdec = otherpartialdec.thisptr
        self.thisptr.copyPartitionStatistics(cpp_otherpartialdec)

    def deleteEmptyBlocks(PartialDecomposition self, bool variables):
        """Deletes empty blocks and sets nblocks accordingly

        A block is considered to be empty if no constraint is assigned to it,
        variables in blocks with no constraints become open

        :param variables: if True, then blocks with no constraints but at least one variable are considered to be nonempty.
        """
        cdef bool cpp_variables = variables
        self.thisptr.deleteEmptyBlocks(cpp_variables)

    def deleteOpencons(PartialDecomposition self, int opencons):
        """Deletes a cons from list of open conss

        :param opencons: id of the cons that is not considered open anymore.
        """
        cdef int cpp_opencons = opencons
        self.thisptr.deleteOpencons(cpp_opencons)

    def deleteOpenvar(PartialDecomposition self, int openvar):
        """Deletes a var from the list of open vars

        :param openvar: id of the var that is not considered open anymore.
        """
        cdef int cpp_openvar = openvar
        self.thisptr.deleteOpenvar(cpp_openvar)

    def displayInfo(PartialDecomposition self, int detailLevel):
        """displays the relevant information of the partialdec

        :param detailLevel: pass a value that indicates how detailed the output should be:
        0: brief overview
        1: block and detector info
        2: cons and var assignments.
        """
        cdef int cpp_detailLevel = detailLevel
        self.thisptr.displayInfo(cpp_detailLevel)

    # def filloutBorderFromConstoblock(PartialDecomposition self, SCIP_HASHMAP * constoblock, int givenNBlocks):
    #     """every constraint is either assigned to master or open

    #     Assignment happens according to the cons assignment information given in constoblock hashmap,
    #     variables are set accordingly
    #     :note: precondition: no constraint or variable is already assigned to a block
    #     :return: scip return code.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def filloutPartialdecFromConstoblock(PartialDecomposition self, SCIP_HASHMAP * constoblock, int givenNBlocks):
    #     """assigns all conss to master or a block

    #     Assignment happens according to the cons assignment information given in constoblock hashmap

    #     :return: scip return code
    #     calculates implicit variable assignment through cons assignment
    #     :note: precondition: no cons or var is already assigned to a block and constoblock contains information for every cons.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def findVarsLinkingToMaster(PartialDecomposition self):
        """reassigns linking variables to master if appropriate

        Variables are reassigned as master if the variable only hits master conss.
        """
        self.thisptr.findVarsLinkingToMaster()

    def findVarsLinkingToStairlinking(PartialDecomposition self):
        """reassigns variables classified as linking to stairlinking if appropriate

        Variables are reassigned as master if the variable hits conss in exactly two consecutive
        blocks.
        """
        self.thisptr.findVarsLinkingToStairlinking()

    def getAncestorID(PartialDecomposition self, int ancestorindex):
        """gets partialdec id of given ancestor id

        :return: partialdec id of given ancestor id.
        """
        cdef int cpp_ancestorindex = ancestorindex
        cdef int result = self.thisptr.getAncestorID(cpp_ancestorindex)
        return result

    def getAncestorList(PartialDecomposition self):
        """get ancestor ids as vector

        :return: vector of ids of all ancestors id.
        """
        cdef vector[int] result = self.thisptr.getAncestorList()
        return result

    def setAncestorList(PartialDecomposition self, object newlist):
        """set ancestor list directly

        :param newlist: new list of ancestor ids.
        """
        cdef vector[int] cpp_newlist = newlist
        self.thisptr.setAncestorList(cpp_newlist)

    def removeAncestorID(PartialDecomposition self, int ancestorid):
        """removes ancestor id from list.
        """
        cdef int cpp_ancestorid = ancestorid
        self.thisptr.removeAncestorID(cpp_ancestorid)

    def addAncestorID(PartialDecomposition self, int ancestor):
        """adds ancestor id to back of list

        :param ancestor: id of ancestor that is to be added.
        """
        cdef int cpp_ancestor = ancestor
        self.thisptr.addAncestorID(cpp_ancestor)

    def getBlocksForRep(PartialDecomposition self, int repid):
        """get a vector of block ids that are identical to block with id repid

        :param repid: id of the representative block
        :return: vector of block ids that are identical to block with id repid.
        """
        cdef int cpp_repid = repid
        cdef vector[int] result = self.thisptr.getBlocksForRep(cpp_repid)
        return result

    def getDetectorClockTime(PartialDecomposition self, int detectorchainindex):
        """returns the time that the detector related to the given detectorchainindex needed for detecting

        :return: the clock time for the corresponding detector in the chain.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getDetectorClockTime(cpp_detectorchainindex)
        return result

    def getDetectorClockTimes(PartialDecomposition self):
        """returns a vector of the clock times that each detector needed that was involved in this partialdec

        :return: vector of the clock times.
        """
        cdef vector[double] result = self.thisptr.getDetectorClockTimes()
        return result

    def getBlockConss(PartialDecomposition self, int block):
        """returns array containing constraints assigned to a block

        :param block: id of the block the constraint indices are returned
        :return: array containing constraints assigned to a block.
        .. note:: This calls the corresponding `getConssForBlock()` method. For consistancy, the method was renamed in the Python API.
        """
        cdef int cpp_block = block
        cdef vector[int] result = self.thisptr.getConssForBlock(cpp_block)
        return [self.getDetprobdata().getCons(c) for c in result]

    # def getDetectorchain(PartialDecomposition self):
    #     """returns detector chain as vector of detector pointers

    #     :return: detector chain as array of detector pointers.
    #     """
    #     cdef vector[DEC_DETECTOR *] result = self.thisptr.getDetectorchain()
    #     return result

    def getFinishedByFinisher(PartialDecomposition self):
        """returns True iff this partialdec was finished by finishPartialdec() method of a detector

        :return: True iff this partialdec was finished by finishPartialdec() method of a detector.
        """
        cdef bool result = self.thisptr.getFinishedByFinisher()
        return result

    def getHashValue(PartialDecomposition self):
        """returns the calculated hash value of this partialdec

        :return: the calculated hash value of this partialdec.
        """
        cdef unsigned long result = self.thisptr.getHashValue()
        return result

    def getID(PartialDecomposition self):
        """returns the unique id of the partialdec

        :return: the unique id of the partialdec.
        """
        cdef int result = self.thisptr.getID()
        return result

    def getLinkingvars(PartialDecomposition self):
        """returns array containing all linking vars indices

        :return: vector containing all linking vars indices
        .. note:: when accessed it is supposed to be sorted.
        """
        cdef vector[int] result = self.thisptr.getLinkingvars()
        return result

    def getMasterconss(PartialDecomposition self):
        """Gets array containing all master conss indices

        :return: array containing all master conss indices
        .. note:: when accessed it is supposed to be sorted.
        """
        cdef vector[int] result = self.thisptr.getMasterconss()
        return [self.getDetprobdata().getCons(c) for c in result]

    def getMastervars(PartialDecomposition self):
        """Gets array containing all master vars indices

        master vars hit only constraints in the master, aka static variables
        :return: array containing all master vars indices.
        """
        cdef vector[int] result = self.thisptr.getMastervars()
        return result

    def getNCoeffsForBlock(PartialDecomposition self, int blockid):
        """Gets the number of nonzero coeffs in a certain block

        :param blockid: of the block the number of nozerors are requested for
        :return: number of nonzero coeffs in a certain block.
        """
        cdef int cpp_blockid = blockid
        cdef int result = self.thisptr.getNCoeffsForBlock(cpp_blockid)
        return result

    def getNCoeffsForMaster(PartialDecomposition self):
        """!Gets the number of nonzero coeffs in master
        :return: the number of nonzero coeffs in master.
        """
        cdef int result = self.thisptr.getNCoeffsForMaster()
        return result

    # def getScore(PartialDecomposition self, SCORETYPE type):
    #     """returns the score of the partialdec (depending on used scoretype)

    #     :param type: the scoretype
    #     :return: the score
    #     @see enum scoretype in cons_decomp.

    #     h
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def hasSetppccardMaster(PartialDecomposition self):
        """checks if all master constraints set partitioning, set packing, set cover, or cardinality constraints

        :return: True iff all master constraints set partitioning, set packing, set cover, or cardinality constraints.
        """
        cdef unsigned int result = self.thisptr.hasSetppccardMaster()
        return result

    def hasSetppcMaster(PartialDecomposition self):
        """checks iff all master constraints set partitioning, set packing, or set cover constraints

        :return: True iff all master constraints set partitioning, set packing, or set cover.
        """
        cdef unsigned int result = self.thisptr.hasSetppcMaster()
        return result

    def hasSetppMaster(PartialDecomposition self):
        """checks iff all master constraints set partitioning, or set packing constraints

        :return: True iff all master constraints set partitioning, or set packing constraints.
        """
        cdef unsigned int result = self.thisptr.hasSetppMaster()
        return result

    def getUsergiven(PartialDecomposition self):
        """Gets the USERGIVEN status of this partialdecs

        :return: the USERGIVEN status of this partialdecs
        """
        # TODO implement function
        raise NotImplementedError()

    def getNAncestors(PartialDecomposition self):
        """Gets number of ancestor partialdecs

        :return: number of ancestor partialdecs.
        """
        cdef int result = self.thisptr.getNAncestors()
        return result

    def getNBlocks(PartialDecomposition self):
        """Gets the number of blocks

        :return: number of blocks.
        """
        cdef int result = self.thisptr.getNBlocks()
        return result

    def getNConss(PartialDecomposition self):
        """Gets the number of constraints

        :return: number of constraints.
        """
        cdef int result = self.thisptr.getNConss()
        return result

    def getNConssForBlock(PartialDecomposition self, int block):
        """Gets size of the vector containing conss assigned to a block

        :param block: id of the block the number of constraints is asked for
        :return: size of the vector containing conss assigned to a block.
        """
        cdef int cpp_block = block
        cdef int result = self.thisptr.getNConssForBlock(cpp_block)
        return result

    def getDetectorchainInfo(PartialDecomposition self):
        """Gets the detectorchain info vector

        :return: detectorchain info vector.
        """
        cdef vector[string] result = self.thisptr.getDetectorchainInfo()
        return result

    def getNDetectors(PartialDecomposition self):
        """Gets the number of detectors the partialdec is propagated by

        :return: number of detectors the partialdec is propagated by.
        """
        cdef int result = self.thisptr.getNDetectors()
        return result

    def getNLinkingvars(PartialDecomposition self):
        """Gets size of the vector containing linking vars

        :return: size of the vector containing linking vars.
        """
        cdef int result = self.thisptr.getNLinkingvars()
        return result

    def getNMasterconss(PartialDecomposition self):
        """Gets size of the vector containing master conss

        :return: size of the vector containing master conss.
        """
        cdef int result = self.thisptr.getNMasterconss()
        return result

    def getNMastervars(PartialDecomposition self):
        """Gets size of the vector containing master vars

        master vars hit only constraints in the master
        :return: size of the vector containing master vars.
        """
        cdef int result = self.thisptr.getNMastervars()
        return result

    def getNNewBlocks(PartialDecomposition self, int detectorchainindex):
        """Gets number of blocks a detector added

        :return: number of blocks a detector added.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef int result = self.thisptr.getNNewBlocks(cpp_detectorchainindex)
        return result

    def getNNewBlocksVector(PartialDecomposition self):
        """gets number of blocks the detectors in the detectorchain added

        :return: number of blocks the detectors in the detectorchain added.
        """
        cdef vector[int] result = self.thisptr.getNNewBlocksVector()
        return result

    def getNTotalStairlinkingvars(PartialDecomposition self):
        """Gets total number of stairlinking vars

        :return: total number of stairlinking vars.
        """
        cdef int result = self.thisptr.getNTotalStairlinkingvars()
        return result

    def getNOpenconss(PartialDecomposition self):
        """Gets size of vector containing constraints not assigned yet

        :return: returns size of vector containing constraints not assigned yet.
        """
        cdef int result = self.thisptr.getNOpenconss()
        return result

    def getNOpenvars(PartialDecomposition self):
        """Gets size of vector containing variables not assigned yet

        :return: size of vector containing variables not assigned yet.
        """
        cdef int result = self.thisptr.getNOpenvars()
        return result

    def getNReps(PartialDecomposition self):
        """Gets the number of blockrepresentatives

        :return: the number of blockrepresentatives.
        """
        cdef int result = self.thisptr.getNReps()
        return result

    def getNStairlinkingvars(PartialDecomposition self, int block):
        """Gets size of the vector containing stairlinking vars

        :param block: id of the block the size of the stairlinking vector is asked for
        :return: size of the vector containing stairlinking vars.
        """
        cdef int cpp_block = block
        cdef int result = self.thisptr.getNStairlinkingvars(cpp_block)
        return result

    def getNVars(PartialDecomposition self):
        """Gets number of vars

        :return: number of vars.
        """
        cdef int result = self.thisptr.getNVars()
        return result

    def getNVarsForBlock(PartialDecomposition self, int block):
        """Gets size of the vector containing vars assigned to a block

        :param block: id of the block the number of variables is asked for
        :return: size of the vector containing vars assigned to a block.
        """
        cdef int cpp_block = block
        cdef int result = self.thisptr.getNVarsForBlock(cpp_block)
        return result

    def getNVarsForBlocks(PartialDecomposition self):
        """Gets overall number of vars assigned to a block

        :return: number of vars that are assigned to any block.
        """
        cdef int result = self.thisptr.getNVarsForBlocks()
        return result

    # def getOpenconss(PartialDecomposition self):
    #     """Gets array containing constraints not assigned yet

    #     :return: array containing constraints not assigned yet.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def getOpenconssVec(PartialDecomposition self):
        """Gets a vector containing constraint ids not assigned yet as vector

        :return: returns a vector containing constraint ids not assigned yet as vector.
        """
        cdef vector[int] result = self.thisptr.getOpenconssVec()
        return result

    def getOpenvars(PartialDecomposition self):
        """Gets array containing variables not assigned yet

        :return: returns array containing variables not assigned yet.
        """
        # TODO implement function
        raise NotImplementedError()

    def getOpenvarsVec(PartialDecomposition self):
        """!Gets array containing variables not assigned yet as vector

        :return: array containing variables not assigned yet as vector.
        """
        cdef vector[int] result = self.thisptr.getOpenvarsVec()
        return result

    def getPctVarsToBorder(PartialDecomposition self, int detectorchainindex):
        """Gets fraction of variables assigned to the border for a detector

        :return: fraction of variables assigned to the border for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsToBorder(cpp_detectorchainindex)
        return result

    def getPctVarsToBorderVector(PartialDecomposition self):
        """Gets fraction of variables assigned to the border for detectors in detectorchain

        :return: vector of fractions of variables assigned to the border for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctVarsToBorderVector()
        return result

    def getPctVarsToBlock(PartialDecomposition self, int detectorchainindex):
        """Gets fraction of variables assigned to a block for a detector

        :return: fraction of variables assigned to a block for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsToBlock(cpp_detectorchainindex)
        return result

    def getPctVarsToBlockVector(PartialDecomposition self):
        """returns fraction of variables assigned to a block for detectors in detectorchain

        :return: vector of fractions of variables assigned to a block for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctVarsToBlockVector()
        return result

    def getPctVarsFromFree(PartialDecomposition self, int detectorchainindex):
        """Gets fraction of variables that are not longer open for a detector

        :return: index of the detector in the detectorchain.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsFromFree(cpp_detectorchainindex)
        return result

    def getPctVarsFromFreeVector(PartialDecomposition self):
        """Gets fraction of variables that are not longer open for detectors in detectorchain

        :return: vector or fractions of variables that are not longer open for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctVarsFromFreeVector()
        return result

    def getPctConssToBorder(PartialDecomposition self, int detectorchainindex):
        """Gets fraction of constraints assigned to the border for a detector

        :return: returns fraction of constraints assigned to the border for a detector
        /
        /**.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssToBorder(cpp_detectorchainindex)
        return result

    def getPctConssToBorderVector(PartialDecomposition self):
        """Gets fraction of constraints assigned to the border for detectors in detectorchain

        :return: vector of fractions of constraints assigned to the border for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctConssToBorderVector()
        return result

    def getPctConssToBlock(PartialDecomposition self, int detectorchainindex):
        """Gets fraction of constraints assigned to a block for a detector

        :return: fraction of constraints assigned to a block for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssToBlock(cpp_detectorchainindex)
        return result

    def getPctConssToBlockVector(PartialDecomposition self):
        """Gets fraction of constraints assigned to a block for detectors in detectorchain

        :return: vector of fractions of constraints assigned to a block for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctConssToBlockVector()
        return result

    def getPctConssFromFree(PartialDecomposition self, int detectorchainindex):
        """Gets fraction of constraints that are not longer open for a detector

        :return: fraction of constraints that are not longer open for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssFromFree(cpp_detectorchainindex)
        return result

    def getPctConssFromFreeVector(PartialDecomposition self):
        """Gets fraction of constraints that are not longer open for detectors in detectorchain

        :return: vector of fractions of constraints that are not longer open for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctConssFromFreeVector()
        return result

    def getRepForBlock(PartialDecomposition self, int blockid):
        """Gets index of the representative block for a block, this might be blockid itself

        :param blockid: id of the block the representative is asked for
        :return: index of the representative block for a block, this might be blockid itself.
        """
        cdef int cpp_blockid = blockid
        cdef int result = self.thisptr.getRepForBlock(cpp_blockid)
        return result

    def getRepVarmap(PartialDecomposition self, int repid, int blockrepid):
        """Gets the represenation varmap

        Var map is vector for represenative repid and the blockrepid-th block that is represented by repid

        :param repid: id of representative
        :param blockrepid: id of block
        :return: the represenation varmap as vector for represenative repid and the blockrepid-th block that is represented by repid.
        """
        cdef int cpp_repid = repid
        cdef int cpp_blockrepid = blockrepid
        cdef vector[int] result = self.thisptr.getRepVarmap(cpp_repid, cpp_blockrepid)
        return result

    def getDetprobdata(PartialDecomposition self):
        """Gets the corresponding detprobdata

        :return: corresponding detprobdata.
        """
        cdef DETPROBDATA * result = self.thisptr.getDetprobdata()
        return DetProbData.create(result)

    def getStairlinkingvars(PartialDecomposition self, int block):
        """Gets array containing stairlinking vars,

        .. note:: if a stairlinking variable links block i and i+1 it is only stored in vector of block i
        :param block: id of the block the stairlinking variable varctor is asked for
        :return: array containing stairlinking vars,.
        """
        cdef int cpp_block = block
        # TODO implement function
        raise NotImplementedError()

    def getVarsForBlock(PartialDecomposition self, int block):
        """Gets array containing vars of a block

        :param block: id of the block the vars are requested for
        :return: returns array containing vars of a block.
        """
        cdef int cpp_block = block
        cdef vector[int] result = self.thisptr.getVarsForBlock(cpp_block)
        return result

    def getVarProbindexForBlock(PartialDecomposition self, int varid, int block):
        """ Gets index in variables array of a block for a variable

        :param varid: the id of the variable the index
        :param block: the corresponding block id
        :return:  returns index in variables array of a block for a variable.
        """
        cdef int cpp_varid = varid
        cdef int cpp_block = block
        cdef int result = self.thisptr.getVarProbindexForBlock(cpp_varid, cpp_block)
        return result

    def isComplete(PartialDecomposition self):
        """Gets whether this partialdec is complete, i.e. it has no more open constraints and variables

        :return: True iff this partialdec is complete
        """
        cdef bool result = self.thisptr.isComplete()
        return result

    def isConsMastercons(PartialDecomposition self, int cons):
        """Gets whether the cons is a master cons

        :param cons: id of ccons to check if it is master constraint
        :return: True iff the cons is a master cons.
        """
        cdef int cpp_cons = cons
        cdef bool result = self.thisptr.isConsMastercons(cpp_cons)
        return result

    def isConsOpencons(PartialDecomposition self, int cons):
        """Gets whether the cons is an open cons

        :param cons: id of cons to check
        :return: True iff the cons is an open cons.
        """
        cdef int cpp_cons = cons
        cdef bool result = self.thisptr.isConsOpencons(cpp_cons)
        return result

    def isAssignedToOrigProb(PartialDecomposition self):
        """Gets whether the partialdec is from the presolved problem

        :return: True iff the partialdec is from the presolved problem.
        """
        cdef bool result = self.thisptr.isAssignedToOrigProb()
        return result

    @property
    def isSelected(PartialDecomposition self):
        """!Gets whether the partialdec is currently selected in explore menue

        :return: True iff the partialdec is currently selected in explore menue.
        """
        cdef bool result = self.thisptr.isSelected()
        return result

    @isSelected.setter
    def isSelected(PartialDecomposition self, bool selected):
        """set the selection status of this partialdecs

        :param selected: whether the partialdec is selected.
        """
        cdef bool cpp_selected = selected
        self.thisptr.setSelected(cpp_selected)

    # def isEqual(PartialDecomposition self, PartialDecomposition otherpartialdec, unsigned int * isequal, bool sortpartialdecs):
    #     """method to check whether this partialdec is equal to a given other partialdec ( \see  isEqual(PARTIALDECOMP*))

    #     :return: scip return code.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def isPropagatedBy(PartialDecomposition self, DEC_DETECTOR * detector):
    #     """Gets whether this partialdec was propagated by specified detector

    #     :param detector: pointer to detector to check for
    #     :return: True iff this partialdec was propagated by detectorID.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def isTrivial(PartialDecomposition self):
        """Gets whether this partialdec is considered to be trivial

        PARTIALDECOMP is considered trivial if all conss are in one block, all conss are in border,
        all variables linking or mastervars, or all constraints and variables are open
        :return: True iff this partialdec is considered to be trivial.
        """
        cdef bool result = self.thisptr.isTrivial()
        return result

    def isVarBlockvarOfBlock(PartialDecomposition self, int var, int block):
        """Checks whether the var is assigned to the block

        :param var: id of var to check
        :param block: id of block to check
        :return: True iff the var is assigned to the block.
        """
        cdef int cpp_var = var
        cdef int cpp_block = block
        cdef bool result = self.thisptr.isVarBlockvarOfBlock(cpp_var, cpp_block)
        return result

    def isVarLinkingvar(PartialDecomposition self, int var):
        """Checks whether the var is a linking var

        :param var: id of var to check
        :return: True iff the var is a linking var.
        """
        cdef int cpp_var = var
        cdef bool result = self.thisptr.isVarLinkingvar(cpp_var)
        return result

    def isVarMastervar(PartialDecomposition self, int var):
        """Checks whether the var is a master var

        :param var: id of var to check
        :return: True iff the var is a master var.
        """
        cdef int cpp_var = var
        cdef bool result = self.thisptr.isVarMastervar(cpp_var)
        return result

    def isVarOpenvar(PartialDecomposition self, int var):
        """Checks whether the var is an open var

        :param var: id of var to check
        :return: True iff the var is an open var
        /
        /**.
        """
        cdef int cpp_var = var
        cdef bool result = self.thisptr.isVarOpenvar(cpp_var)
        return result

    def isVarStairlinkingvar(PartialDecomposition self, int var):
        """Checks whether the var is a stairlinking var

        :param var: id of var to check
        :return: True iff the var is a stairlinking var.
        """
        cdef int cpp_var = var
        cdef bool result = self.thisptr.isVarStairlinkingvar(cpp_var)
        return result

    def isVarStairlinkingvarOfBlock(PartialDecomposition self, int var, int block):
        """Checks whether the var is a stairlinkingvar of a specified block

        :param var: id of var to check if it is a stairlinking variable hitting specified block
        :param block: id of block to check
        :return: True iff the var is a stairlinkingvar of a specified block.
        """
        cdef int cpp_var = var
        cdef int cpp_block = block
        cdef bool result = self.thisptr.isVarStairlinkingvarOfBlock(cpp_var, cpp_block)
        return result

    # def printPartitionInformation(PartialDecomposition self, SCIP * givenscip, FILE * file):
    #     """prints partition information as described in \see cls reader

    #     :param givenscip: scip data structure
    #     :param file: output file.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def refineToBlocks(PartialDecomposition self):
        """refine partialdec with focus on blocks

        strategy: assigns open conss and vars if they can be found in blocks
        (without respect to open vars and conss  see assignHittingOpenconss(), see assignHittingOpenvars())
        .. note:: partialdec might be not complete.
        """
        self.thisptr.refineToBlocks()

    def refineToMaster(PartialDecomposition self):
        """refine partialdec with focus on master

        strategy: do obvious ( see considerImplicits()) assignments and
        assign other conss and vars to master if possible (see assignOpenPartialHittingToMaster()).
        """
        self.thisptr.refineToMaster()

    def setConsPartitionStatistics(PartialDecomposition self, int detectorchainindex, ConsPart partition, object consclassesmaster):
        """registers statistics for a used conspartition.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef ConsPartition* cpp_partition = partition.consPartition
        cdef vector[int] cpp_consclassesmaster = consclassesmaster
        self.thisptr.setConsPartitionStatistics(cpp_detectorchainindex, cpp_partition, cpp_consclassesmaster)

    def setConsToBlock(PartialDecomposition self, int consToBlock, int block):
        """adds a constraint to a block, does not delete this cons from list of open conss

        :param consToBlock: id of cons to add
        :param block: id of block to add.
        """
        cdef int cpp_consToBlock = consToBlock
        cdef int cpp_block = block
        self.thisptr.setConsToBlock(cpp_consToBlock, cpp_block)

    # def fixConsToBlock(PartialDecomposition self, int cons, int block):
    #     """adds a constraint to a block

    #     :param cons: id of cons to add
    #     :param block: id of block to add.
    #     """
    #     cdef int cpp_cons = cons
    #     cdef int cpp_block = block
    #     self.thisptr.fixConsToBlock(cpp_cons, cpp_block)

    def setConsToMaster(PartialDecomposition self, int consToMaster):
        """adds a constraint to the master constraints, does not delete this cons from list of open conss

        :param consToMaster: id of cons to add.
        """
        cdef int cpp_consToMaster = consToMaster
        self.thisptr.setConsToMaster(cpp_consToMaster)

    # def fixConsToMaster(PartialDecomposition self, int cons):
    #     """fixes a constraint to the master constraints

    #     :param cons: id of cons to add
    #     @warning This method modifies the vector PARTIALDECOMP::openconss! Hence, any kind of iterator might be invalid afterwards!.
    #     """
    #     cdef int cpp_cons = cons
    #     self.thisptr.fixConsToMaster(cpp_cons)

    # def setDetectorchain(PartialDecomposition self, object givenDetectorChain):
    #     """sets the detectorchain with the given vector of detector pointers

    #     :param givenDetectorChain: vector of detector pointers.
    #     """
    #     cdef vector[DEC_DETECTOR *] cpp_givenDetectorChain = givenDetectorChain
    #     self.thisptr.setDetectorchain(cpp_givenDetectorChain)

    # def setDetectorPropagated(PartialDecomposition self, DEC_DETECTOR * detector):
    #     """sets partialdec to be propagated by a detector

    #     :param detector: pointer to detector that is registered for this partialdec.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def setDetectorFinished(PartialDecomposition self, DEC_DETECTOR * detector):
    #     """sets detector that finished the partialdec

    #     :param detector: pointer to detector that has finished this partialdecs.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def setDetectorFinishedOrig(PartialDecomposition self, DEC_DETECTOR * detectorID):
    #     """sets detector that finished the partialdec in the original problem

    #     :param detectorID: pointer to detector that has finished this partialdecs
    #     :note: does not add the detector to the detectorchain and does not modify partition statistics.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def setFinishedByFinisher(PartialDecomposition self, bool finished):
        """sets whether this partialdec was finished by a finishing detector

        :param finished: is this partialdecs finished by a finishing detector.
        """
        cdef bool cpp_finished = finished
        self.thisptr.setFinishedByFinisher(cpp_finished)

    def setFinishedByFinisherOrig(PartialDecomposition self, bool finished):
        """sets whether this partialdec was finished by a finishing detector in the original problem (in case this partialdec was translated)

        :param finished: was this partialdecs finished by a finishing detector in orig.
        """
        cdef bool cpp_finished = finished
        self.thisptr.setFinishedByFinisherOrig(cpp_finished)

    def setNBlocks(PartialDecomposition self, int nblocks):
        """sets number of blocks, only increasing number allowed

        :param nblocks: new number of blocks.
        """
        cdef int cpp_nblocks = nblocks
        self.thisptr.setNBlocks(cpp_nblocks)

    def setStemsFromOrig(PartialDecomposition self, bool fromorig):
        """sets whether this partialdec stems from an orig problem partialdec

        :param fromorig: has this partialdec ancestors from the orig problem.
        """
        cdef bool cpp_fromorig = fromorig
        self.thisptr.setStemsFromOrig(cpp_fromorig)

    # def setUsergiven(PartialDecomposition self, cpp.USERGIVEN usergiven):
    #     """sets whether this partialdec is user given

    #     :param usergiven: is this partialdec user given.
    #     """
    #     self.thisptr.setUsergiven(usergiven)

    def setVarPartitionStatistics(PartialDecomposition self, int detectorchainindex, VarPart partition, object varclasseslinking, object varclassesmaster):
        """registers statistics for a used varpartition.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef VarPartition * cpp_partition = partition.thisptr
        cdef vector[int] cpp_varclasseslinking = varclasseslinking
        cdef vector[int] cpp_varclassesmaster = varclassesmaster
        self.thisptr.setVarPartitionStatistics(cpp_detectorchainindex, cpp_partition, cpp_varclasseslinking, cpp_varclassesmaster)

    def setVarToBlock(PartialDecomposition self, int varToBlock, int block):
        """adds a variable to the linking variables, does not delete this var from list of open vars

        :param varToBlock: id of var to be added
        :param block: id of block to be added.
        """
        cdef int cpp_varToBlock = varToBlock
        cdef int cpp_block = block
        self.thisptr.setVarToBlock(cpp_varToBlock, cpp_block)

    def fixVarToBlock(PartialDecomposition self, int var, int block):
        """adds a variable to the linking variables

        :param var: id of var to be added
        :param block: id of block to be added.
        """
        cdef int cpp_var = var
        cdef int cpp_block = block
        self.thisptr.fixVarToBlock(cpp_var, cpp_block)

    def setVarToLinking(PartialDecomposition self, int varToLinking):
        """adds a variable to the linking variables, does not delete this var from list of open vars

        :param varToLinking: var to be set to linking.
        """
        cdef int cpp_varToLinking = varToLinking
        self.thisptr.setVarToLinking(cpp_varToLinking)

    def fixVarToLinking(PartialDecomposition self, int var):
        """adds a variable to the linking variables

        :param var: var to be set to linking.
        """
        cdef int cpp_var = var
        self.thisptr.fixVarToLinking(cpp_var)

    def setVarToMaster(PartialDecomposition self, int varToMaster):
        """adds a variable to the master variables, does not delete this var from list of open vars

        master variables hit only constraints in the master.
        """
        cdef int cpp_varToMaster = varToMaster
        self.thisptr.setVarToMaster(cpp_varToMaster)

    def fixVarToMaster(PartialDecomposition self, int var):
        """adds a variable to the master variables

        master variables hit only constraints in the master.
        """
        cdef int cpp_var = var
        self.thisptr.fixVarToMaster(cpp_var)

    def setVarToStairlinking(PartialDecomposition self, int varToStairLinking, int block1, int block2):
        """adds a variable to the stairlinking variables, does not delete this var from list of open vars

        :param varToStairLinking: id of variable to be added
        :param block1: id of block one
        :param block2: id of block two
        .. note:: stairlinking variables are only registered in block with smaller index.
        """
        cdef int cpp_varToStairLinking = varToStairLinking
        cdef int cpp_block1 = block1
        cdef int cpp_block2 = block2
        self.thisptr.setVarToStairlinking(cpp_varToStairLinking, cpp_block1, cpp_block2)

    def fixVarToStairlinking(PartialDecomposition self, int var, int firstblock):
        """adds a variable to the stairlinking variables

        :param var: id of variable to be added
        :param firstblock: stairlinking variables hit exactly two consecutive blocks, this is the index of the first of these blocks
        .. note:: stairlinking variables are only registered in block with smaller index.
        """
        cdef int cpp_var = var
        cdef int cpp_firstblock = firstblock
        self.thisptr.fixVarToStairlinking(cpp_var, cpp_firstblock)

    def fixConsToBlockByName(PartialDecomposition self, consname, int blockid):
        """assigns a constraint by name to a block

        .. seealso:: * :meth:`fixConsToBlock`
        :return: True iff successful.
        """
        c_consname = str_conversion(consname)
        cdef int cpp_blockid = blockid
        cdef bool result = self.thisptr.fixConsToBlockByName(c_consname, cpp_blockid)
        return result

    def fixVarToBlockByName(PartialDecomposition self, varname, int blockid):
        """assigns a variable by name to a block

        .. seealso:: * :meth:`fixVarToBlock`
        :return: True iff successful.
        """
        c_varname = str_conversion(varname)
        cdef int cpp_blockid = blockid
        cdef bool result = self.thisptr.fixVarToBlockByName(c_varname, cpp_blockid)
        return result

    def fixConsToMasterByName(PartialDecomposition self, consname):
        """assgins a constraint by name as master

        .. seealso:: * :meth:`fixConsToMaster`
        :return: True iff successful.
        """
        c_consname = str_conversion(consname)
        cdef bool result = self.thisptr.fixConsToMasterByName(c_consname)
        return result

    def fixVarToMasterByName(PartialDecomposition self, varname):
        """assigns a variable with given name as master

        .. seealso:: * :meth:`fixVarToMaster`
        :return: True iff successful.
        """
        c_varname = str_conversion(varname)
        cdef bool result = self.thisptr.fixVarToMasterByName(c_varname)
        return result

    def fixVarToLinkingByName(PartialDecomposition self, varname):
        """assigns a variable by name to the linking variables

        .. seealso:: * :meth:`fixVarToLinking`
        :return: True iff successful.
        """
        c_varname = str_conversion(varname)
        cdef bool result = self.thisptr.fixVarToLinkingByName(c_varname)
        return result

    def showVisualization(PartialDecomposition self):
        """generates and opens a gp visualization of the partialdec

        .. note:: linux only.
        """
        self.thisptr.showVisualization()

    # def generateVisualization(PartialDecomposition self, filename, outname, GP_OUTPUT_FORMAT outputformat):
    #     """generates a visualization of the partialdec using gnuplot

    #     :param filename: Path where to store the gp file
    #     :param outname: Path at which gnuplot will output its result
    #     :param outputformat: The format of the gnuplot output.

    #     Should match the file extension of outname
    #     :note: linux only, requires gnuplot
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def writeVisualizationFile(PartialDecomposition self, filename, outname, GP_OUTPUT_FORMAT outputformat):
    #     """writes a gp visualization of the partialdec to a file

    #     :param filename: Path where to store the gp file
    #     :param outname: Path at which gnuplot will output its result
    #     :param outputformat: The format of the gnuplot output.

    #     Should match the file extension of outname
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def shouldCompletedByConsToMaster(PartialDecomposition self):
        """Checks whether this partialdec is a userpartialdec that should be completed

        the completion should be done by setting unspecified constraints to master
        :return: True iff this partialdec is a userpartialdec that should be completed.
        """
        cdef unsigned int result = self.thisptr.shouldCompletedByConsToMaster()
        return result

    def sort(PartialDecomposition self):
        """sorts the vars and conss data structures by their indices

        :return: True if the internal order of variables or constraints changed.
        """
        cdef bool result = self.thisptr.sort()
        return result

    def setPctConssToBlockVector(PartialDecomposition self, object newvector):
        """set statistical vector of fractions of constraints set to blocks per involved detector

        :param newvector: vector of fractions of constraints set to blocks per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssToBlockVector(cpp_newvector)

    def setPctConssFromFreeVector(PartialDecomposition self, object newvector):
        """set statistical vector of fractions of constraints that are not longer open per involved detector

        :param newvector: vector of fractions of constraints that are not longer open per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssFromFreeVector(cpp_newvector)

    def setPctConssToBorderVector(PartialDecomposition self, object newvector):
        """set statistical vector of fractions of constraints assigned to the border per involved detector

        :param newvector: vector of fractions of constraints assigned to the border per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssToBorderVector(cpp_newvector)

    def setPctVarsToBorderVector(PartialDecomposition self, object newvector):
        """set statistical vector of fraction of variables assigned to the border per involved detector

        :param newvector: vector of fractions of variables assigned to the border per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsToBorderVector(cpp_newvector)

    def setPctVarsToBlockVector(PartialDecomposition self, object newvector):
        """set statistical vector of fractions of variables assigned to a block per involved detector

        :param newvector: vector of fractions of variables assigned to a block per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsToBlockVector(cpp_newvector)

    def setPctVarsFromFreeVector(PartialDecomposition self, object newvector):
        """set statistical vector of variables that are not longer open per involved detector

        :param newvector: vector of fractions of variables that are not longer open per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsFromFreeVector(cpp_newvector)

    def setDetectorClockTimes(PartialDecomposition self, object newvector):
        """set statistical vector of the times that the detectors needed for detecting per involved detector

        :param newvector: vector of the times that the detectors needed for detecting per involved detector.
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setDetectorClockTimes(cpp_newvector)

    @property
    def classicScore(PartialDecomposition self):
        """gets the classic score

        .. note:: -1 iff not calculated yet, .. seealso:: GCGconshdlrDecompCalcClassicScore
        :return: border area score.
        """
        cdef double result = self.thisptr.getClassicScore()
        return result

    @classicScore.setter
    def classicScore(PartialDecomposition self, double score):
        """set the classic score.
        """
        cdef double cpp_score = score
        self.thisptr.setClassicScore(cpp_score)

    @property
    def borderAreaScore(PartialDecomposition self):
        """gets the border area score

        .. note:: -1 iff not calculated yet, .. seealso:: GCGconshdlrDecompCalcBorderAreaScore
        :return: border area score.
        """
        cdef double result = self.thisptr.getBorderAreaScore()
        return result

    @borderAreaScore.setter
    def borderAreaScore(PartialDecomposition self, double score):
        """set the border area score.
        """
        cdef double cpp_score = score
        self.thisptr.setBorderAreaScore(cpp_score)

    @property
    def maxWhiteScore(PartialDecomposition self):
        """gets the maximum white area score

        "maximum white score" is fraction of the area of the decomposed matrix that is neither block or border
        .. note:: -1 iff not calculated yet, .. seealso:: GCGconshdlrDecompCalcMaxWhiteScore
        :return: maximum  white area score
        """
        cdef double result = self.thisptr.getMaxWhiteScore()
        return result

    @maxWhiteScore.setter
    def maxWhiteScore(PartialDecomposition self, double score):
        """set the maximum white area score.
        """
        cdef double cpp_score = score
        self.thisptr.setMaxWhiteScore(cpp_score)

    @property
    def maxForWhiteScore(PartialDecomposition self):
        """gets the maximum foreseeing white area score

        .. note:: -1 iff not calculated yet, .. seealso:: GCGconshdlrDecompCalcMaxForseeingWhiteScore
        :return: maximum foreseeing white area score
        """
        cdef double result = self.thisptr.getMaxForWhiteScore()
        return result

    @maxForWhiteScore.setter
    def maxForWhiteScore(PartialDecomposition self, double score):
        """set the maximum foreseeing white area score.
        """
        cdef double cpp_score = score
        self.thisptr.setMaxForWhiteScore(cpp_score)

    @property
    def partForWhiteScore(PartialDecomposition self):
        """gets the setpartitioning maximum foreseeing white area score

        .. note:: -1 iff not calculated yet, .. seealso:: GGCGconshdlrDecompCalcSetPartForseeingWhiteScore
        :return: setpartitioning maximum foreseeing white area score
        """
        cdef double result = self.thisptr.getSetPartForWhiteScore()
        return result

    @partForWhiteScore.setter
    def partForWhiteScore(PartialDecomposition self, double score):
        """set the setpartitioning maximum foreseeing white area score.
        """
        cdef double cpp_score = score
        self.thisptr.setSetPartForWhiteScore(cpp_score)

    @property
    def maxForWhiteAggScore(PartialDecomposition self):
        """gets the maximum foreseeing white area score with respect to aggregatable blocks

        .. note:: -1 iff not calculated yet, .. seealso:: GCGconshdlrDecompCalcMaxForeseeingWhiteAggScore
        :return: maximum foreseeing white area score with respect to aggregatable blocks
        """
        cdef double result = self.thisptr.getMaxForWhiteAggScore()
        return result

    @maxForWhiteAggScore.setter
    def maxForWhiteAggScore(PartialDecomposition self, double score):
        """set the maximum foreseeing white area score with respect to aggregatable blocks.
        """
        cdef double cpp_score = score
        self.thisptr.setMaxForWhiteAggScore(cpp_score)

    @property
    def partForWhiteAggScore(PartialDecomposition self):
        """gets the setpartitioning maximum foreseeing white area score with respect to aggregateable

        .. note:: -1 iff not calculated yet, .. seealso:: GCGconshdlrDecompCalcSetPartForWhiteAggScore
        :return: setpartitioning maximum foreseeing white area score with respect to aggregateable.
        """
        cdef double result = self.thisptr.getSetPartForWhiteAggScore()
        return result

    @partForWhiteAggScore.setter
    def partForWhiteAggScore(PartialDecomposition self, double score):
        """set the setpartitioning maximum foreseeing white area score with respect to aggregateable.
        """
        cdef double cpp_score = score
        self.thisptr.setSetPartForWhiteAggScore(cpp_score)

    @property
    def bendersScore(PartialDecomposition self):
        """gets the benders score

        .. note:: -1 iff not calculated yet, .. seealso:: GCGconshdlrDecompCalcBendersScore
        :return: benders score.
        """
        cdef double result = self.thisptr.getBendersScore()
        return result

    @bendersScore.setter
    def bendersScore(PartialDecomposition self, double score):
        """set the benders score.
        """
        cdef double cpp_score = score
        self.thisptr.setBendersScore(cpp_score)

    @property
    def strongDecompScore(PartialDecomposition self):
        """gets the strong decomposition score

        .. note:: -1 iff not calculated yet, .. seealso:: GCGconshdlrDecompCalcStrongDecompositionScore
        :return: strong decomposition score.
        """
        cdef double result = self.thisptr.getStrongDecompScore()
        return result

    @strongDecompScore.setter
    def strongDecompScore(PartialDecomposition self, double score):
        """set the strong decomposition score.
        """
        cdef double cpp_score = score
        self.thisptr.setStrongDecompScore(cpp_score)

    def prepare(PartialDecomposition self):
        """sorts the partialdec and calculates a its implicit assignments, hashvalue and evaluation

        :return: SCIP_OKAY if the result is consistent, SCIP_ERROR if there was an inconsistency.
        """
        self.thisptr.prepare()

    def aggInfoCalculated(PartialDecomposition self):
        """Checks if the aggregation information was already calculated

        :return: True iff the aggregation information was already calculated.
        """
        cdef bool result = self.thisptr.aggInfoCalculated()
        return result

    def calcAggregationInformation(PartialDecomposition self, bool ignoreDetectionLimits):
        """computes if aggregation of sub problems is possible

        checks if aggregation of sub problems is possible and stores the corresponding aggregation information

        :param ignoreDetectionLimits: Set to True if computation should ignore detection limits.

        This parameter is ignored if the patched bliss version is not present.
        """
        cdef bool cpp_ignoreDetectionLimits = ignoreDetectionLimits
        self.thisptr.calcAggregationInformation(cpp_ignoreDetectionLimits)

    def getConssForBlocks(PartialDecomposition self):
        cdef vector[vector[int]] result = self.thisptr.getConssForBlocks()
        return result

    def getTranslatedpartialdecid(PartialDecomposition self):
        cdef int result = self.thisptr.getTranslatedpartialdecid()
        return result

    def setTranslatedpartialdecid(PartialDecomposition self, int decid):
        cdef int cpp_decid = decid
        self.thisptr.setTranslatedpartialdecid(cpp_decid)

    def buildDecChainString(PartialDecomposition self, buffer):
        """creates a detector chain short string for this partialdec, is built from detector chain.
        """
        c_buffer = str_conversion(buffer)
        self.thisptr.buildDecChainString(c_buffer)

    # END AUTOGENERATED BLOCK

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
                else:
                    c_output_format = GP_OUTPUT_FORMAT_SVG #as default output format

                self.thisptr.generateVisualization(c_gp_filename, c_outfile, c_output_format)

                if format == "svg":
                    data = outfile.read_text()
                elif format == "png":
                    data = outfile.read_bytes()
                self._visualizations[format] = data

        return self._visualizations[format]

    def __repr__(PartialDecomposition self):
        return f"<PartialDecomposition: nBlocks={self.getNBlocks()}, nMasterConss={self.getNMasterconss()}, nMasterVars={self.getNMastervars()}, nLinkingVars={self.getNLinkingvars()}, maxForWhiteScore={self.maxForWhiteScore}>"
