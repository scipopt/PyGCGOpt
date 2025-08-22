cdef class PartialDecomposition:
    """class to manage partial decompositions

    each partialdec corresponds to one :class:`DetProbData` which contains the problem information,
    there is one detprobdata for the original and the transformed problem.
    """
    cdef PARTIALDECOMP* thisptr
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

    def __copy__(self):
        return PartialDecomposition.create(new PARTIALDECOMP(self.thisptr))

    def copy(self):
        return copy(self)

    @property
    def isSelected(self):
        """gets whether the partialdec is currently selected in explore menue

        :return: True iff the partialdec is currently selected in explore menue.
        """
        cdef bool result = self.thisptr.isSelected()
        return result

    @isSelected.setter
    def isSelected(self, bool selected):
        """set the selection status of this partialdecs

        :param selected: whether the partialdec is selected
        """
        cdef bool cpp_selected = selected
        self.thisptr.setSelected(cpp_selected)

    def addBlock(self):
        """adds a new block

        :return: the number (id) of the new block
        """
        cdef int result = self.thisptr.addBlock()
        return result

    def assignOpenConssToMaster(self):
        """assigns open constraints to master
        """
        self.thisptr.assignOpenConssToMaster()

    def isComplete(self):
        """gets whether this partialdec is complete, i.e. it has no more open constraints and variables

        :return: True iff this partialdec is complete
        """
        cdef bool result = self.thisptr.isComplete()
        return result

    def isTrivial(self):
        """gets whether this partialdec is considered to be trivial

        :return: True iff this partialdec is considered to be trivial

        .. note:: PartialDecomposition is considered trivial if all constraints are in one block, all constraints are in border,
        all variables linking or mastervariables, or all constraints and variables are open
        """
        cdef bool result = self.thisptr.isTrivial()
        return result

    def getNBlocks(self):
        """gets the number of blocks

        :return: number of blocks
        """
        cdef int result = self.thisptr.getNBlocks()
        return result

    def getNConss(self):
        """gets the number of constraints

        :return: number of constraints
        """
        cdef int result = self.thisptr.getNConss()
        return result

    def getNVars(self):
        """gets number of variables

        :return: number of variables
        """
        cdef int result = self.thisptr.getNVars()
        return result

    def getID(self):
        """gets the unique id of the partialdecomposition

        :return: unique id of the partialdecomposition
        """
        cdef int result = self.thisptr.getID()
        return result

    def getHashValue(self):
        """gets the calculated hash value of the partialdecomposition

        :return: calculated hash value of the partialdecomposition
        """
        cdef unsigned long result = self.thisptr.getHashValue()
        return result

    def getDetprobdata(self):
        """gets the corresponding detprobdata

        :return: corresponding detprobdata
        """
        cdef DETPROBDATA* result = self.thisptr.getDetprobdata()
        return DetProbData.create(result)

    def fixConsToMaster(self, Constraint cons):
        """fixes a Constraint to the master constraints

        :param cons: scip#Constraint to add
        :type cons: Constraint
        """
        self.thisptr.fixConsToMaster(cons.scip_cons)

    def fixConssToMaster(self, conss):
        """fixes all constraints to the master constraints.

        :param conss: An iterable of scip#Constraint objects
        :type conss: An iterable of scip#Constraint objects
        :raises TypeError: occurs if conss is not an Iterable
        .. seealso:: * :meth:`fixConsToMaster`
        """
        if not isinstance(conss, Iterable):
            raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(conss)))
        for cons in conss:
            self.fixConsToMaster(<Constraint?>cons)

    def fixConsToBlockId(self, Constraint cons, int block_id):
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
        :param block_id: id of block to add
        .. seealso:: * :meth:`fixConsToBlockId`
                     * :meth:`fixConssToBlock`
        """
        if not isinstance(conss, Iterable):
            raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(conss)))
        for cons in conss:
            self.fixConsToBlockId(<Constraint?>cons, block_id)

    #def fixVarToBlockId(self, Variable var, int block_id):
    #    """adds a variable to the block

    #    :param var: variable to be added
    #    :param block_id: id of block to be added
    #    """
    #    self.thisptr.fixVarToBlock(var.scip_var, block_id)

    #def fixVarsToBlockId(self, vars, int block_id):
    #    """adds all variables to the block

    #    :param vars: An iterable of scip#Variable objects
    #    :param block_id: id of block to be added
    #    """
    #    if not isinstance(vars, Iterable):
    #        raise TypeError("Expected iterable as first argument. Got '{}' instead.".format(type(vars)))
    #    for var in vars:
    #        self.fixVarToBlockId(<Variable?>var, block_id)

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

    def hasSetppccardMaster(self):
        """checks if all master constraints set partitioning, set packing, set cover or cardinality constraints

        :return: True iff all master constraints set partitioning, set packing, set cover or cardinality constraints
        """
        cdef bool result = self.thisptr.hasSetppccardMaster()
        return result

    def hasSetppcMaster(self):
        """checks if all master constraints set partitioning, set packing or set cover constraints

        :return: True iff all master constraints set partitioning, set packing or set cover
        """
        cdef bool result = self.thisptr.hasSetppcMaster()
        return result

    def hasSetppMaster(self):
        """checks if all master constraints set partitioning or set packing constraints

        :return: True iff all master constraints set partitioning or set packing constraints
        """
        cdef bool result = self.thisptr.hasSetppMaster()
        return result

    def getNOpenconss(self):
        """gets number of constraints that are opened (not assigned yet)

        :return: number of open constraints
        """
        cdef int result = self.thisptr.getNOpenconss()
        return result

    def getOpenconss(self):
        """gets a list of constraints that are opened (not assigned yet)

        :return: list of constraints
        """
        cdef vector[int] result = self.thisptr.getOpenconssVec()
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        return [Constraint.create(det_prob_data.getCons(consIndex)) for consIndex in result]

    def getNLinkingvars(self):
        """gets number of linking variables

        :return: number of linking variables
        """
        cdef int result = self.thisptr.getNLinkingvars()
        return result

    def getLinkingvars(self):
        """gets a list of variables that are assigned to linking variables

        :return: list of variables
        """
        cdef vector[int] result = self.thisptr.getLinkingvars()
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        return [Variable.create(det_prob_data.getVar(varIndex)) for varIndex in result]

    def getNMasterconss(self):
        """gets number of master constraints

        :return: number of master constraints
        """
        cdef int result = self.thisptr.getNMasterconss()
        return result

    def getMasterconss(self):
        """gets a list of constraints that are assigned to master variables

        :return: list of constraints
        """
        cdef vector[int] result = self.thisptr.getMasterconss()
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        return [Constraint.create(det_prob_data.getCons(consIndex)) for consIndex in result]

    def getNMastervars(self):
        """gets number of master variables

        :return: number of master variables
        """
        cdef int result = self.thisptr.getNMastervars()
        return result

    def getMastervars(self):
        """gets a list of variables that are assigned to master variables (static variables)

        :return: list of variables
        """
        cdef vector[int] result = self.thisptr.getMastervars()
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        return [Variable.create(det_prob_data.getVar(varIndex)) for varIndex in result]

    def getNOpenvars(self):
        """gets number of variables that are opened (not assigned yet)

        :return: number of open variables
        """
        cdef int result = self.thisptr.getNOpenvars()
        return result

    def getOpenvars(self):
        """gets a list of variables that are opened (not assigned yet)

        :return: list of variables
        """
        cdef vector[int] result = self.thisptr.getOpenvarsVec()
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        return [Variable.create(det_prob_data.getVar(varIndex)) for varIndex in result]

    def getNStairlinkingvars(self, int block_id):
        """gets number of variables that are assigned to stairlinking variables of the specified block

        :param block_id: id of the block the number of stairlinking variables is asked for
        :return: number of stairlinking variables

        .. note:: if a stairlinking variable links block i and i+1 it is only stored in vector of block i
        """
        cdef int result = self.thisptr.getNStairlinkingvars(block_id)
        return result

    #def getStairlinkingvars(self, int block_id):
        """gets a list of variables that are assigned to stairlinking variables of the specified block

        :param block_id: id of the block the stairlinking variables is asked for
        :return: list of variables

        .. note:: if a stairlinking variable links block i and i+1 it is only stored in vector of block i
        """
    #    cdef vector[int] result = self.thisptr.getStairlinkingvarsVec(block_id)
    #    cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
    #    return [Variable.create(det_prob_data.getVar(varIndex)) for varIndex in result]

    def getConssForBlock(self, int block_id):
        """gets a list of constraints that are assigned to the specified block

        :param block_id: id of the block the constraints is asked for
        :return: list of constraints
        """
        cdef vector[int] result = self.thisptr.getConssForBlock(block_id)
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        return [Constraint.create(det_prob_data.getCons(consIndex)) for consIndex in result]

    def getConssForBlocks(self):
        cdef vector[vector[int]] result = self.thisptr.getConssForBlocks()
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        return [[Constraint.create(det_prob_data.getCons(consIndex)) for consIndex in block] for block in result]

    def getNVarsForBlock(self, int block_id):
        """gets number of variables that are assigned to the specified block

        :param block_id: id of the block the number of variables is asked for
        :return: number of variables
        """
        cdef int result = self.thisptr.getNVarsForBlock(block_id)
        return result

    def getVarsForBlock(self, int block_id):
        """gets a list of variables that are assigned to the specified block

        :param block_id: id of the block the variables is asked for
        :return: list of variables
        """
        cdef vector[int] result = self.thisptr.getVarsForBlock(block_id)
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        return [Variable.create(det_prob_data.getVar(varIndex)) for varIndex in result]

    def setConsToMaster(self, Constraint cons):
        """adds constraint to the master constraints, does not delete this constraint from list of open constraints

        :param cons: scip#Constraint to add
        :type cons: Constraint
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.setConsToMaster(det_prob_data.getIndexForCons(cons.scip_cons))

    def setVarToBlockId(self, Variable var, int block_id):
        """adds variable to the linking variables, does not delete this variable from list of open variables

        :param var: scip#Variable to add
        :type var: Variable
        :param block_id: id of block to add
        :type block_id: int
        """
        cdef int cpp_block_id = block_id
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.setVarToBlock(det_prob_data.getIndexForVar(var.scip_var), cpp_block_id)

    def setVarToLinking(self, Variable var):
        """adds variable to the linking variables, does not delete this variable from list of open variables

        :param var: scip#Variable to add
        :type var: Variable
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.setVarToLinking(det_prob_data.getIndexForVar(var.scip_var))

    def fixVarToLinking(self, Variable var):
        """adds variable to the linking variables

        :param var: scip#Variable to add
        :type var: Variable
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.fixVarToLinking(det_prob_data.getIndexForVar(var.scip_var))

    def setVarToMaster(self, Variable var):
        """adds variable to the master variables, does not delete this variable from list of open variables

        master variables hit only constraints in the master

        :param var: scip#Variable to add
        :type var: Variable
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.setVarToMaster(det_prob_data.getIndexForVar(var.scip_var))

    def fixVarToMaster(self, Variable var):
        """adds variable to the master variables

        master variables hit only constraints in the master

        :param var: scip#Variable to add
        :type var: Variable
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.fixVarToMaster(det_prob_data.getIndexForVar(var.scip_var))

    def setVarToStairlinking(self, Variable var, int block1_id, int block2_id):
        """adds variable to the stairlinking variables, does not delete this variable from list of open variables

        :param var: scip#Variable to add
        :type var: Variable
        :param block1_id: id of block one
        :type block1_id: int
        :param block2_id: id of block two
        :type block2_id: int
        .. note:: stairlinking variables are only registered in block with smaller index
        """
        cdef int cpp_block1_id = block1_id
        cdef int cpp_block2_id = block2_id
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.setVarToStairlinking(det_prob_data.getIndexForVar(var.scip_var), cpp_block1_id, cpp_block2_id)

    def fixVarToStairlinking(self, Variable var, int firstblock_id):
        """adds variable to the stairlinking variables

        :param var: scip#Variable to add
        :type var: Variable
        :param firstblock_id: stairlinking variables hit exactly two consecutive blocks, this is the index of the first of these blocks
        :type firstblock_id: int
        .. note:: stairlinking variables are only registered in block with smaller index
        """
        cdef int cpp_firstblock_id = firstblock_id
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.fixVarToStairlinking(det_prob_data.getIndexForVar(var.scip_var), cpp_firstblock_id)

    def isConsMastercons(self, Constraint cons):
        """gets whether the constraint is a master constraint

        :param cons: constraint to check if it is master constraint
        :type cons: scip#Constraint
        :return: True iff the constraint is a master constraint
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        cdef bool result = self.thisptr.isConsMastercons(det_prob_data.getIndexForCons(cons.scip_cons))
        return result

    def isConsOpencons(self, Constraint cons):
        """gets whether the constraint is an open constraint

        :param cons: constraint to check if it is open constraint
        :type cons: scip#Constraint
        :return: True iff the constraint is an open constraint
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        cdef bool result = self.thisptr.isConsOpencons(det_prob_data.getIndexForCons(cons.scip_cons))
        return result

    def isVarBlockvarOfBlock(self, Variable var, int block_id):
        """gets whether the variable is assigned to the block

        :param var: variable to check if it is in the specified block
        :type var: scip#Variable
        :param block_id: id of block to check
        :return: True iff the variable is assigned to the specified block
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        cdef bool result = self.thisptr.isVarBlockvarOfBlock(det_prob_data.getIndexForVar(var.scip_var), block_id)
        return result

    def isVarLinkingvar(self, Variable var):
        """gets whether the variable is a linking variable

        :param var: variable to check if it is in a linking variable
        :type var: scip#Variable
        :return: True iff the variable is a linking var
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        cdef bool result = self.thisptr.isVarLinkingvar(det_prob_data.getIndexForVar(var.scip_var))
        return result

    def isVarMastervar(self, Variable var):
        """gets whether the variable is a master variable

        :param var: variable to check if it is a master variable
        :type var: scip#Variable
        :return: True iff the variable is a master variable
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        cdef bool result = self.thisptr.isVarMastervar(det_prob_data.getIndexForVar(var.scip_var))
        return result

    def isVarOpenvar(self, Variable var):
        """gets whether the variable is an open variable

        :param var: variable to check if it is an open variable
        :type var: scip#Variable
        :return: True iff the variable is an open variable
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        cdef bool result = self.thisptr.isVarOpenvar(det_prob_data.getIndexForVar(var.scip_var))
        return result

    def isVarStairlinkingvar(self, Variable var):
        """gets whether the variable is a stairlinking variable

        :param var: variable to check if it is a stairlinking variable
        :type var: scip#Variable
        :return: True iff the variable is a stairlinking variable
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        cdef bool result = self.thisptr.isVarStairlinkingvar(det_prob_data.getIndexForVar(var.scip_var))
        return result

    def isVarStairlinkingvarOfBlock(self, Variable var, int block_id):
        """checks whether the var is a stairlinkingvar of a specified block

        :param var: variable to check if it is a stairlinking variable hitting the specified block
        :type var: scip#Variable
        :param block_id: id of block to check
        :return: True iff the variable is a stairlinking variable of the specified block
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        cdef bool result = self.thisptr.isVarStairlinkingvarOfBlock(det_prob_data.getIndexForVar(var.scip_var), block_id)
        return result

    def deleteOpencons(self, Constraint cons):
        """deletes constraint from list of open constraints

        :param cons: scip#Constraint that is not considered open anymore
        :type cons: Constraint
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.deleteOpencons(det_prob_data.getIndexForCons(cons.scip_cons))

    def deleteOpenvar(self, Variable var):
        """deletes variable from the list of open variables

        :param var: scip#Variable that is not considered open anymore
        :type var: Variable
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.deleteOpenvar(det_prob_data.getIndexForVar(var.scip_var))

    def setUsergiven(self, USERGIVEN value=COMPLETED_CONSTOMASTER):
        self.thisptr.setUsergiven(value)

    def getUsergiven(self):
        """gets the PY_USERGIVEN status of this partialdecs

        :return: the PY_USERGIVEN status of this partialdecs
        """
        return self.thisptr.getUsergiven()

    @property
    def max_white_score(self):
        """gets the maximum white area score

        "maximum white score" is fraction of the area of the decomposed matrix that is neither block or border
        .. note:: -1 iff not calculated yet

        :return: maximum  white area score
        """
        return self.thisptr.getMaxWhiteScore()

    @property
    def classic_score(self):
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion("classic")))

    @property
    def border_area_score(self):
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion("border area")))

    @property
    def max_for_white_score(self):
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion("max foreseeing white")))

    @property
    def set_part_for_white_score(self):
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion("ppc-max-white with aggregation info")))

    @property
    def max_for_white_agg_score(self):
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion("max foreseeing white with aggregation info")))

    @property
    def set_part_for_white_agg_score(self):
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion("ppc-max-white")))

    @property
    def benders_score(self):
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion("experimental benders score")))

    @property
    def strong_decomp_score(self):
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion("strong decomposition score")))

    def addClockTime(self, double clocktime):
        """adds detection time of one detector

        incorporates the needed time of some detector in the detector chain.
        """
        cdef double cpp_clocktime = clocktime
        self.thisptr.addClockTime(cpp_clocktime)

    def addPctConssFromFree(self, double pct):
        """adds percentage of closed constraints

        bookkeeping information: fraction of constraints that are not longer open for a detector added to detector chain
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssFromFree(cpp_pct)

    def addPctConssToBlock(self, double pct):
        """adds percentage of constraints assigned to blocks

        bookkeeping information: adds fraction of constraints assigned to a block for a detector added to detector chain
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssToBlock(cpp_pct)

    def addPctConssToBorder(self, double pct):
        """adds percentage of constraints assigned to border

        bookkeeping information: adds fraction of constraints assigned to the border for a detector added to detector chain
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctConssToBorder(cpp_pct)

    def addPctVarsFromFree(self, double pct):
        """adds percentage of closed variables

        bookkeeping information: adds fraction of variables that are not longer open for a detector added to detector chain
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsFromFree(cpp_pct)

    def addPctVarsToBlock(self, double pct):
        """adds percentage of variables assigned to blocks

        bookkeeping information: adds fraction of variables assigned to a block for a detector added to detector chain
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsToBlock(cpp_pct)

    def addPctVarsToBorder(self, double pct):
        """adds percentage of variables assigned to border

        bookkeeping information: adds fraction of variables assigned to the border for a detector added to detector chain
        """
        cdef double cpp_pct = pct
        self.thisptr.addPctVarsToBorder(cpp_pct)

    def getPctVarsToBorder(self, int detectorchainindex):
        """gets fraction of variables assigned to the border for a detector

        :return: fraction of variables assigned to the border for a detector
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsToBorder(cpp_detectorchainindex)
        return result

    def getPctVarsToBorderVector(self):
        """gets fraction of variables assigned to the border for detectors in detectorchain

        :return: vector of fractions of variables assigned to the border for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctVarsToBorderVector()
        return result

    def getPctVarsToBlock(self, int detectorchainindex):
        """gets fraction of variables assigned to a block for a detector

        :return: fraction of variables assigned to a block for a detector
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsToBlock(cpp_detectorchainindex)
        return result

    def getPctVarsToBlockVector(self):
        """returns fraction of variables assigned to a block for detectors in detectorchain

        :return: vector of fractions of variables assigned to a block for detectors in detectorchain
        """
        cdef vector[double] result = self.thisptr.getPctVarsToBlockVector()
        return result

    def getPctVarsFromFree(self, int detectorchainindex):
        """Gets fraction of variables that are not longer open for a detector

        :return: index of the detector in the detectorchain
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctVarsFromFree(cpp_detectorchainindex)
        return result

    def getPctVarsFromFreeVector(self):
        """gets fraction of variables that are not longer open for detectors in detectorchain

        :return: vector or fractions of variables that are not longer open for detectors in detectorchain
        """
        cdef vector[double] result = self.thisptr.getPctVarsFromFreeVector()
        return result

    def getPctConssToBorder(self, int detectorchainindex):
        """gets fraction of constraints assigned to the border for a detector

        :return: returns fraction of constraints assigned to the border for a detector
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssToBorder(cpp_detectorchainindex)
        return result

    def getPctConssToBorderVector(self):
        """gets fraction of constraints assigned to the border for detectors in detectorchain

        :return: vector of fractions of constraints assigned to the border for detectors in detectorchain
        """
        cdef vector[double] result = self.thisptr.getPctConssToBorderVector()
        return result

    def getPctConssToBlock(self, int detectorchainindex):
        """gets fraction of constraints assigned to a block for a detector

        :return: fraction of constraints assigned to a block for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssToBlock(cpp_detectorchainindex)
        return result

    def getPctConssToBlockVector(self):
        """gets fraction of constraints assigned to a block for detectors in detectorchain

        :return: vector of fractions of constraints assigned to a block for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctConssToBlockVector()
        return result

    def getPctConssFromFree(self, int detectorchainindex):
        """gets fraction of constraints that are not longer open for a detector

        :return: fraction of constraints that are not longer open for a detector.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef double result = self.thisptr.getPctConssFromFree(cpp_detectorchainindex)
        return result

    def getPctConssFromFreeVector(self):
        """gets fraction of constraints that are not longer open for detectors in detectorchain

        :return: vector of fractions of constraints that are not longer open for detectors in detectorchain.
        """
        cdef vector[double] result = self.thisptr.getPctConssFromFreeVector()
        return result

    def getScore(self, Score score):
        """gets value of the score

        :param score: score for which the value should be returned
        :return: value of score
        """
        return self.thisptr.getScore(GCGfindScore(self.thisptr.getDetprobdata().getScip(), str_conversion(score.scorename)))

    def calcStairlinkingVars(self):
        """reassigns linking variables to stairlinkingvars if possible

        potentially reorders blocks for making a maximum number of linking vars stairlinking
        if all vars that connect exactly two blocks have a staircase structure, all of them become stairlinkingvars
        otherwise, the stairlinking assignment is done greedily
        .. note:: precondition: partialdec does not have any stairlinking vars.
        """
        self.thisptr.calcStairlinkingVars()

    def checkAllConssAssigned(self):
        """checks if all constraints are assigned and deletes the open constraint vector if so

        :return: True iff all constraints are assigned
        """
        cdef bool result = self.thisptr.checkAllConssAssigned()
        return result

    def checkConsistency(self):
        """checks whether the assignments in the partialdec are consistent

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

    def complete(self):
        """assigns all open constraints and open variables trivially

        strategy: assigns all open conss and vars to blocks if they can be refined there, otherwise to the master

        :note: partialdecomps should usually be completed by a detector, only use this function if you know what you are doing
        """
        self.thisptr.complete()

    def sort(self):
        """sorts the variables and constraints data structures by their indices

        :return: True if the internal order of variables or constraints changed
        """
        cdef bool result = self.thisptr.sort()
        return result

    def completeByConnected(self):
        """assigns all open constraints and open variables

        strategy: assigns all constraints and variables to the same block if they are connected,
        a constraint and a variable are adjacent if the variable appears in the constraint
        """
        self.thisptr.completeByConnected()

    def completeByConnectedConssAdjacency(self):
        """assigns all open constraints and open variables

        strategy: assigns all constraints and variables to the same block if they are connected
        a constraint and a variable are adjacent if the variable appears in the constraint
        .. note:: this relies on the consadjacency structure of the detprobdata
        hence it cannot be applied in presence of linking variables.
        """
        self.thisptr.completeByConnectedConssAdjacency()

    def completeGreedily(self):
        """assigns all open constraints and open variables

        strategy: assigns a constraint (and related variables) to a new block if possible,
        if not to an existing block if possible (by means of prior var assignments)
        and finally to master, if there does not exist such a block.
        """
        self.thisptr.completeGreedily()

    def prepare(self):
        """sorts the partialdec and calculates a its implicit assignments, hashvalue and evaluation
        
        """
        self.thisptr.prepare()

    def getNDetectors(self):
        """gets the number of detectors the partialdec is propagated by

        :return: number of detectors the partialdec is propagated by
        """
        cdef int result = self.thisptr.getNDetectors()
        return result

    def getNReps(self):
        """gets the number of blockrepresentatives

        :return: the number of blockrepresentatives
        """
        cdef int result = self.thisptr.getNReps()
        return result

    def getNCoeffsForMaster(self):
        """gets number of nonzero coefficients in master

        :return: number of nonzero coefficients in master
        """
        cdef int result = self.thisptr.getNCoeffsForMaster()
        return result

    def setNBlocks(self, int nblocks):
        """sets number of blocks, only increasing number allowed

        :param nblocks: new number of blocks
        """
        cdef int cpp_nblocks = nblocks
        self.thisptr.setNBlocks(cpp_nblocks)

    def getNTotalStairlinkingvars(self):
        """gets total number of stairlinking variables

        :return: total number of stairlinking variables
        """
        cdef int result = self.thisptr.getNTotalStairlinkingvars()
        return result

    def getNVarsForBlocks(self):
        """gets overall number of variables assigned to a block

        :return: number of variables that are assigned to any block
        """
        cdef int result = self.thisptr.getNVarsForBlocks()
        return result

    def findVarsLinkingToMaster(self):
        """reassigns linking variables to master if appropriate

        Variables are reassigned as master if the variable only hits master conss
        """
        self.thisptr.findVarsLinkingToMaster()

    def findVarsLinkingToStairlinking(self):
        """reassigns variables classified as linking to stairlinking if appropriate

        Variables are reassigned as master if the variable hits conss in exactly two consecutive blocks
        """
        self.thisptr.findVarsLinkingToStairlinking()

    def getNCoeffsForBlock(self, int block_id):
        """gets number of nonzero coeffs in a certain block

        :param block_id: of the block the number of nozerors are requested for
        :return: number of nonzero coeffs in a certain block
        """
        cdef int cpp_block_id = block_id
        cdef int result = self.thisptr.getNCoeffsForBlock(cpp_block_id)
        return result

    def addNNewBlocks(self, int nnewblocks):
        """adds how many new blocks were introduced

        bookkeeping information: adds number of new blocks created by a detector added to detector chain
        """
        cdef int cpp_nnewblocks = nnewblocks
        self.thisptr.addNNewBlocks(cpp_nnewblocks)

    def getNConssForBlock(self, int block_id):
        """gets size of the vector containing constraints assigned to a block

        :param block_id: id of the block the number of constraints is asked for
        :return: size of the vector containing constraints assigned to a block
        """
        cdef int cpp_block_id = block_id
        cdef int result = self.thisptr.getNConssForBlock(cpp_block_id)
        return result

    def removeMastercons(self, Constraint cons):
        """removes the given constraint from master

        :param cons: constraint to be removed from master
        :type cons: scip#Constraint
        """
        cdef DETPROBDATA* det_prob_data = self.thisptr.getDetprobdata()
        self.thisptr.removeMastercons(det_prob_data.getIndexForCons(cons.scip_cons))

    def isAssignedToOrigProb(self):
        """gets whether the partialdec is from the presolved problem

        :return: True iff the partialdec is from the presolved problem
        """
        cdef bool result = self.thisptr.isAssignedToOrigProb()
        return result

    def refineToBlocks(self):
        """refine partialdec with focus on blocks

        strategy: assigns open conss and vars if they can be found in blocks
        (without respect to open vars and conss  see assignHittingOpenconss(), see assignHittingOpenvars())
        .. note:: partialdec might be not complete.
        """
        self.thisptr.refineToBlocks()

    def refineToMaster(self):
        """refine partialdec with focus on master

        strategy: do obvious ( see considerImplicits()) assignments and
        assign other conss and vars to master if possible (see assignOpenPartialHittingToMaster())
        """
        self.thisptr.refineToMaster()

    def alreadyAssignedConssToBlocks(self):
        """checks if at least one constraint is assigned to some block

        :return: True iff at least one constraint is assigned to a block
        """
        cdef bool result = self.thisptr.alreadyAssignedConssToBlocks()
        return result

    def assignCurrentStairlinking(self):
        """assigns open variables to stairlinking if appropriate

        assigns open variables to stairlinking if they can be found in exactly two consecutive blocks

        :return: True iff at least one stairlinkingvariable was assigned
        """
        cdef bool result = self.thisptr.assignCurrentStairlinking()
        return result

    def aggInfoCalculated(self):
        """checks if the aggregation information was already calculated

        :return: True iff the aggregation information was already calculated
        """
        cdef bool result = self.thisptr.aggInfoCalculated()
        return result

    def considerImplicits(self):
        """assigns every open constraint/variable

        Assignments happen as follows:
        - to the respective block if it hits exactly one blockvariable/blockconstraint and no open variables/constraints
        - to master/linking if it hits blockvariables/blockconstraints assigned to different blocks
        - and every constraint to master that hits a master variable
        - and every variable to master if it does not hit any blockconstraint and has no open constraint
        - leave the constraitn/variable open if nothing from the above holds
        """
        self.thisptr.considerImplicits()

    def displayInfo(self, int detailLevel):
        """displays the relevant information of the partialdec

        :param detailLevel: pass a value that indicates how detailed the output should be:
        0: brief overview
        1: block and detector info
        2: cons and var assignments
        """
        cdef int cpp_detailLevel = detailLevel
        self.thisptr.displayInfo(cpp_detailLevel)

    def getAncestorList(self):
        """gets ancestor ids as vector

        :return: list of ids of all ancestors id
        """
        cdef vector[int] result = self.thisptr.getAncestorList()
        return result

    def assignSmallestComponentsButOneConssAdjacency(self):
        """computes components by connectedness of constraints and variables

        computes components corresponding to connectedness of constraints and variables
        and assigns them accordingly (all but one of largest components)

        strategy: assigns all conss same block if they are connected
        two constraints are adjacent if there is a common variable

        .. note:: this relies on the consadjacency structure of the detprobdata
        hence it cannot be applied in presence of linking variables.
        """
        self.thisptr.assignSmallestComponentsButOneConssAdjacency()

    def deleteEmptyBlocks(self, bool variables):
        """deletes empty blocks and sets nblocks accordingly

        A block is considered to be empty if no constraint is assigned to it,
        variables in blocks with no constraints become open

        :param variables: if True, then blocks with no constraints but at least one variable are considered to be nonempty
        """
        cdef bool cpp_variables = variables
        self.thisptr.deleteEmptyBlocks(cpp_variables)

    def getTranslatedpartialdecid(self):
        cdef int result = self.thisptr.getTranslatedpartialdecid()
        return result

    def setTranslatedpartialdecid(self, int dec_id):
        cdef int cpp_dec_id = dec_id
        self.thisptr.setTranslatedpartialdecid(cpp_dec_id)

    def buildDecChainString(self, buffer):
        """creates a detector chain short string for this partialdec, is built from detector chain

        """
        c_buffer = str_conversion(buffer)
        self.thisptr.buildDecChainString(c_buffer)

    def addDecChangesFromAncestor(self, PartialDecomposition ancestor):
        """adds the statistical differences to an ancestor

        incorporates the changes from ancestor partialdec into the statistical data structures
        """
        cdef PARTIALDECOMP* cpp_ancestor = ancestor.thisptr
        self.thisptr.addDecChangesFromAncestor(cpp_ancestor)

    def addDetectorChainInfo(PartialDecomposition self, decinfo):
        """adds information about the detector chain

        adds a detectorchain information string to the corresponding vector
        (that carries information for each detector call)
        """
        c_decinfo = str_conversion(decinfo)
        self.thisptr.addDetectorChainInfo(c_decinfo)

    def getAncestorID(self, int ancestor_id):
        """gets partialdec id of given ancestor id

        :return: partialdec id of given ancestor id
        """
        cdef int cpp_ancestor_id = ancestor_id
        cdef int result = self.thisptr.getAncestorID(cpp_ancestor_id)
        return result

    def copyPartitionStatistics(PartialDecomposition self, PartialDecomposition otherpartialdec):
        """copies the given partialdec's partition statistics

        :param otherpartialdec: partialdec whose partition statistics are to be copied
        """
        cdef PARTIALDECOMP* cpp_otherpartialdec = otherpartialdec.thisptr
        self.thisptr.copyPartitionStatistics(cpp_otherpartialdec)

    def removeAncestorID(PartialDecomposition self, int ancestor_id):
        """removes ancestor id from list
        """
        cdef int cpp_ancestor_id = ancestor_id
        self.thisptr.removeAncestorID(cpp_ancestor_id)

    def addAncestorID(self, int ancestor_id):
        """adds ancestor id to back of list

        :param ancestor: id of ancestor that is to be added
        """
        cdef int cpp_ancestor_id = ancestor_id
        self.thisptr.addAncestorID(cpp_ancestor_id)

    def getBlocksForRep(PartialDecomposition self, int rep_id):
        """get a vector of block ids that are identical to block with id rep_id

        :param rep_id: id of the representative block
        :return: vector of block ids that are identical to block with id rep_id
        """
        cdef int cpp_rep_id = rep_id
        cdef vector[int] result = self.thisptr.getBlocksForRep(cpp_rep_id)
        return result

    def getNNewBlocksVector(self):
        """gets number of blocks the detectors in the detectorchain added

        :return: number of blocks the detectors in the detectorchain added
        """
        cdef vector[int] result = self.thisptr.getNNewBlocksVector()
        return result

    def getNAncestors(self):
        """gets number of ancestor partialdecs

        :return: number of ancestor partialdecs
        """
        cdef int result = self.thisptr.getNAncestors()
        return result

    def getNNewBlocks(self, int detectorchain_id):
        """gets number of blocks a detector added

        :return: number of blocks a detector added
        """
        cdef int cpp_detectorchain_id = detectorchain_id
        cdef int result = self.thisptr.getNNewBlocks(cpp_detectorchain_id)
        return result

    def getFinishedByFinisher(self):
        """returns True iff this partialdec was finished by finishPartialdec() method of a detector

        :return: True iff this partialdec was finished by finishPartialdec() method of a detector
        """
        cdef bool result = self.thisptr.getFinishedByFinisher()
        return result

    def getDetectorchainInfo(self):
        """gets the detectorchain info vector

        :return: detectorchain info vector
        """
        cdef vector[string] result = self.thisptr.getDetectorchainInfo()
        return result

    def shouldCompletedByConsToMaster(self):
        """checks whether this partialdec is a userpartialdec that should be completed

        the completion should be done by setting unspecified constraints to master
        :return: True iff this partialdec is a userpartialdec that should be completed
        """
        cdef unsigned int result = self.thisptr.shouldCompletedByConsToMaster()
        return result

    def setFinishedByFinisher(self, bool finished):
        """sets whether this partialdec was finished by a finishing detector

        :param finished: is this partialdecs finished by a finishing detector
        """
        cdef bool cpp_finished = finished
        self.thisptr.setFinishedByFinisher(cpp_finished)

    def calcAggregationInformation(self, bool ignoreDetectionLimits):
        """computes if aggregation of sub problems is possible

        checks if aggregation of sub problems is possible and stores the corresponding aggregation information

        :param ignoreDetectionLimits: Set to True if computation should ignore detection limits

        This parameter is ignored if the patched bliss version is not present.
        """
        cdef bool cpp_ignoreDetectionLimits = ignoreDetectionLimits
        self.thisptr.calcAggregationInformation(cpp_ignoreDetectionLimits)

    def setDetectorPropagated(self, Detector detector):
        """sets partialdec to be propagated by a detector

        :param detector: detector that is registered for this partialdec
        """
        cdef GCG_DETECTOR* dec = <GCG_DETECTOR*>detector
        self.thisptr.setDetectorPropagated(dec)

    def setDetectorFinished(self, Detector detector):
        """sets detector that finished the partialdec

        :param detector: detector that has finished this partialdecs
        """
        cdef GCG_DETECTOR* dec = <GCG_DETECTOR*>detector
        self.thisptr.setDetectorFinished(dec)

    def setDetectorFinishedOrig(self, Detector detector):
        """sets detector that finished the partialdec in the original problem

        :param detector: detector that has finished this partialdecs
        :note: does not add the detector to the detectorchain and does not modify partition statistics
        """
        cdef GCG_DETECTOR* dec = <GCG_DETECTOR*>detector
        self.thisptr.setDetectorFinishedOrig(dec)

    def setPctConssToBlockVector(self, object newvector):
        """set statistical vector of fractions of constraints set to blocks per involved detector

        :param newvector: vector of fractions of constraints set to blocks per involved detector
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssToBlockVector(cpp_newvector)

    def setPctConssFromFreeVector(self, object newvector):
        """set statistical vector of fractions of constraints that are not longer open per involved detector

        :param newvector: vector of fractions of constraints that are not longer open per involved detector
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssFromFreeVector(cpp_newvector)

    def setPctConssToBorderVector(self, object newvector):
        """set statistical vector of fractions of constraints assigned to the border per involved detector

        :param newvector: vector of fractions of constraints assigned to the border per involved detector
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctConssToBorderVector(cpp_newvector)

    def setPctVarsToBorderVector(self, object newvector):
        """set statistical vector of fraction of variables assigned to the border per involved detector

        :param newvector: vector of fractions of variables assigned to the border per involved detector
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsToBorderVector(cpp_newvector)

    def setPctVarsToBlockVector(self, object newvector):
        """set statistical vector of fractions of variables assigned to a block per involved detector

        :param newvector: vector of fractions of variables assigned to a block per involved detector
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsToBlockVector(cpp_newvector)

    def setPctVarsFromFreeVector(self, object newvector):
        """set statistical vector of variables that are not longer open per involved detector

        :param newvector: vector of fractions of variables that are not longer open per involved detector
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setPctVarsFromFreeVector(cpp_newvector)

    def setDetectorClockTimes(self, object newvector):
        """set statistical vector of the times that the detectors needed for detecting per involved detector

        :param newvector: vector of the times that the detectors needed for detecting per involved detector
        """
        cdef vector[double] cpp_newvector = newvector
        self.thisptr.setDetectorClockTimes(cpp_newvector)

    def setFinishedByFinisherOrig(self, bool finished):
        """sets whether this partialdec was finished by a finishing detector in the original problem (in case this partialdec was translated)

        :param finished: was this partialdecs finished by a finishing detector in orig
        """
        cdef bool cpp_finished = finished
        self.thisptr.setFinishedByFinisherOrig(cpp_finished)

    def setStemsFromOrig(self, bool fromorig):
        """sets whether this partialdec stems from an orig problem partialdec

        :param fromorig: has this partialdec ancestors from the orig problem
        """
        cdef bool cpp_fromorig = fromorig
        self.thisptr.setStemsFromOrig(cpp_fromorig)

    def getDetectorClockTime(self, int detectorchain_id):
        """returns the time that the detector related to the given detectorchainindex needed for detecting

        :return: the clock time for the corresponding detector in the chain
        """
        cdef int cpp_detectorchain_id = detectorchain_id
        cdef double result = self.thisptr.getDetectorClockTime(cpp_detectorchain_id)
        return result

    def getRepForBlock(self, int block_id):
        """gets index of the representative block for a block, this might be block_id itself

        :param block_id: id of the block the representative is asked for
        :return: index of the representative block for a block, this might be block_id itself
        """
        cdef int cpp_block_id = block_id
        cdef int result = self.thisptr.getRepForBlock(cpp_block_id)
        return result

    def showVisualization(self):
        """generates and opens a gp visualization of the partialdec

        .. note:: linux only
        """
        self.thisptr.showVisualization()

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

    def __repr__(self):
        return f"<PartialDecomposition: nBlocks={self.getNBlocks()}, nMasterConss={self.getNMasterconss()}, nMasterVars={self.getNMastervars()}, nLinkingVars={self.getNLinkingvars()}, max_for_white_score={self.max_for_white_score()}>"















    # BEGIN AUTOGENERATED BLOCK
    def setAncestorList(PartialDecomposition self, object newlist):
        """set ancestor list directly

        :param newlist: new list of ancestor ids.
        """
        cdef vector[int] cpp_newlist = newlist
        self.thisptr.setAncestorList(cpp_newlist)

    def getDetectorClockTimes(PartialDecomposition self):
        """returns a vector of the clock times that each detector needed that was involved in this partialdec

        :return: vector of the clock times.
        """
        cdef vector[double] result = self.thisptr.getDetectorClockTimes()
        return result

    # def getDetectorchain(PartialDecomposition self):
    #     """returns detector chain as vector of detector pointers

    #     :return: detector chain as array of detector pointers.
    #     """
    #     cdef vector[GCG_DETECTOR *] result = self.thisptr.getDetectorchain()
    #     return result

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

    def getVarProbindexForBlock(PartialDecomposition self, int varid, int block):
        """Gets index in variables array of a block for a variable

        :param varid: the id of the variable the index
        :param block: the corresponding block id
        :return:  returns index in variables array of a block for a variable.
        """
        cdef int cpp_varid = varid
        cdef int cpp_block = block
        cdef int result = self.thisptr.getVarProbindexForBlock(cpp_varid, cpp_block)
        return result

    # def isEqual(PartialDecomposition self, PartialDecomposition otherpartialdec, unsigned int * isequal, bool sortpartialdecs):
    #     """method to check whether this partialdec is equal to a given other partialdec ( \see  isEqual(PARTIALDECOMP*))

    #     :return: scip return code.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    # def isPropagatedBy(PartialDecomposition self, GCG_DETECTOR * detector):
    #     """Gets whether this partialdec was propagated by specified detector

    #     :param detector: pointer to detector to check for
    #     :return: True iff this partialdec was propagated by detectorID.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    def setConsPartitionStatistics(PartialDecomposition self, int detectorchainindex, ConsPart partition, object consclassesmaster):
        """registers statistics for a used conspartition.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef ConsPartition* cpp_partition = partition.consPartition
        cdef vector[int] cpp_consclassesmaster = consclassesmaster
        self.thisptr.setConsPartitionStatistics(cpp_detectorchainindex, cpp_partition, cpp_consclassesmaster)

    #def setDetectorchain(PartialDecomposition self, object givenDetectorChain):
    #     """sets the detectorchain with the given vector of detector pointers

    #     :param givenDetectorChain: vector of detector pointers.
    #     """
    #     cdef vector[GCG_DETECTOR *] cpp_givenDetectorChain = givenDetectorChain
    #     self.thisptr.setDetectorchain(cpp_givenDetectorChain)

    def setVarPartitionStatistics(PartialDecomposition self, int detectorchainindex, VarPart partition, object varclasseslinking, object varclassesmaster):
        """registers statistics for a used varpartition.
        """
        cdef int cpp_detectorchainindex = detectorchainindex
        cdef VarPartition * cpp_partition = partition.varPartition
        cdef vector[int] cpp_varclasseslinking = varclasseslinking
        cdef vector[int] cpp_varclassesmaster = varclassesmaster
        self.thisptr.setVarPartitionStatistics(cpp_detectorchainindex, cpp_partition, cpp_varclasseslinking, cpp_varclassesmaster)
