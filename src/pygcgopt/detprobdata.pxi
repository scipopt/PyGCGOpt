cdef class DetProbData:
    """Base class holding a pointer to corresponding DETPROBDATA"""

    @staticmethod
    cdef create(DETPROBDATA* thisptr):
        if thisptr == NULL:
            raise Warning("cannot create DetProbData with DETPROBDATA* == NULL")
        new_DetProbData = DetProbData()
        new_DetProbData.thisptr = thisptr
        return new_DetProbData

    def createConsPart(self, name, nclasses, ncons):
        """returns a ConsPart object

        :param name: name of the constraint partition
        :type name: str
        :param nclasses: number of classes
        :type nclasses: int
        :param ncons: number of constraints
        :type ncons: int
        :return: constraint partition object
        :rtype: :class:`ConsPart`
        """
        cdef SCIP* scip = self.thisptr.getScip()
        c_name = str_conversion(name)
        cdef ConsPartition* conspatitionptr = new ConsPartition(scip, c_name, nclasses, ncons)
        return ConsPart.create(conspatitionptr, self)

    def createVarPart(self, name, nclasses, nvars):
        """returns a VarPart object

        :param name: name of the variable partition
        :type name: str
        :param nclasses: number of classes
        :type nclasses: int
        :param nvars: number of variables
        :type nvars: int
        :return: variable partition object
        :rtype: :class:`VarPart`
        """
        cdef SCIP* scip = self.thisptr.getScip()
        c_name = str_conversion(name)
        cdef VarPartition* varpatitionptr = new VarPartition(scip, c_name, nclasses, nvars)
        return VarPart.create(varpatitionptr, self)
        
    def candidatesNBlocks(self):
        """candidate for the number of blocks, second int indicates how often a candidate was added
        """
        cdef vector[pair[int, int]] result = self.thisptr.candidatesNBlocks
        return result

    def candidatesNBlocks(self, object candidatesNBlocks):
        cdef vector[pair[int, int]] cpp_candidatesNBlocks = candidatesNBlocks
        self.thisptr.candidatesNBlocks = cpp_candidatesNBlocks

    def conspartitioncollection(self):
        """collection of different constraint class distributions.
        """
        cdef vector[ConsPartition *] result = self.thisptr.conspartitioncollection
        return [ConsPart.create(r, <DetProbData>weakref.proxy(self)) for r in result]

    def conspartitioncollection(self, object conspartitioncollection):
        # this seems to be possible only when we use C++11 (-std=c++11)
        # maybe it will be fixed in a future version of Cython
        cdef vector[ConsPartition *] cpp_conspartitioncollection
        cdef ConsPartition * conspartitioncollection_ptr = NULL
        cdef ConsPart conspartitioncollection_element
        for conspartitioncollection_element in conspartitioncollection:
            conspartitioncollection_ptr = <ConsPartition*> conspartitioncollection_element.consPartition
            cpp_conspartitioncollection.push_back(conspartitioncollection_ptr)
        self.thisptr.conspartitioncollection = cpp_conspartitioncollection

    def varpartitioncollection(self):
        """collection of different variable class distributions.
        """
        cdef vector[VarPartition *] result = self.thisptr.varpartitioncollection
        return [VarPart.create(r, self) for r in result]

    def varpartitioncollection(self, object varpartitioncollection):
        # this seems to be possible only when we use C++11 (-std=c++11)
        # maybe it will be fixed in a future version of Cython
        cdef vector[VarPartition *] cpp_varpartitioncollection
        cdef VarPartition * varpartitioncollection_ptr = NULL
        cdef VarPart varpartitioncollection_element
        for varpartitioncollection_element in varpartitioncollection:
            varpartitioncollection_ptr = <VarPartition*> varpartitioncollection_element.varPartition
            cpp_varpartitioncollection.push_back(varpartitioncollection_ptr)
        self.thisptr.varpartitioncollection = cpp_varpartitioncollection

    def classificationtime(self):
        """time that was consumed by the classification of the constraint and variables classifiers.
        """
        cdef double result = self.thisptr.classificationtime
        return result

    def classificationtime(self, double classificationtime):
        cdef double cpp_classificationtime = classificationtime
        self.thisptr.classificationtime = cpp_classificationtime

    def nblockscandidatescalctime(self):
        """time that was used to calulate the candidates of te block number.
        """
        cdef double result = self.thisptr.nblockscandidatescalctime
        return result

    def nblockscandidatescalctime(DetProbData self, double nblockscandidatescalctime):
        cdef double cpp_nblockscandidatescalctime = nblockscandidatescalctime
        self.thisptr.nblockscandidatescalctime = cpp_nblockscandidatescalctime

    def postprocessingtime(self):
        """time that was spent in postproceesing decomposigtions.
        """
        cdef double result = self.thisptr.postprocessingtime
        return result

    def postprocessingtime(self, double postprocessingtime):
        cdef double cpp_postprocessingtime = postprocessingtime
        self.thisptr.postprocessingtime = cpp_postprocessingtime

    def translatingtime(self):
        """time that was spent by transforming partialdecs between presolved and orig problem.
        """
        cdef double result = self.thisptr.translatingtime
        return result

    def translatingtime(self, double translatingtime):
        cdef double cpp_translatingtime = translatingtime
        self.thisptr.translatingtime = cpp_translatingtime

    def addConsPartition(self, ConsPart partition):
        """adds a constraint partition if it is no duplicate of an existing constraint partition

        :param partition: constraint partition object
        :type partition: :class:`ConsPart`
        """
        cdef ConsPartition* cpp_partition = partition.consPartition
        self.thisptr.addConsPartition(cpp_partition)

    def addCandidatesNBlocksNVotes(self, int candidate, int nvotes):
        """adds a candidate for block number and counts how often a candidate is added

        :param candidate: candidate for block size
        :type candidate: int
        :param nvotes: number of votes this candidates will get
        :type nvotes: int
        """
        self.thisptr.addCandidatesNBlocksNVotes(candidate, nvotes)

    def addPartialdecToAncestor(self, PartialDecomposition partialdec):
        """adds a partialdec to ancestor partialdecs

        :param partialdec: partialdec that is added to the ancestor partialdecs
        :type partialdec: :class:`PartialDecomposition`
        """
        cdef PARTIALDECOMP* cpp_partialdec = partialdec.thisptr
        self.thisptr.addPartialdecToAncestor(cpp_partialdec)

    def addPartialdecToOpen(self, PartialDecomposition partialdec):
        """adds a partialdec to current partialdecs (data structure for partialdecs that are going to processed in the propagation rounds)

        :param partialdec: pointer of partialdec to be added
        :type partialdec: :class:`PartialDecomposition`
        :return: True if the partialdecs was successfully added (i.e. it is no duplicate of a known partialdec)
        :rtype: bool
        """
        cdef PARTIALDECOMP* cpp_partialdec = partialdec.thisptr
        return self.thisptr.addPartialdecToOpen(cpp_partialdec)

    def addPartialdecToFinished(self, PartialDecomposition partialdec):
        """adds a partialdec to finished partialdecs

        .. seealso:: * :meth:`addPartialdecToFinishedUnchecked()`

        :param partialdec: pointer of partialdec that is going to be added to the finished partialdecs (data structure to carry finished decompositions)
        :type partialdec: :class:`PartialDecomposition`
        :return: True if the partialdecs was successfully added (i.e. it is no duplicate of a known partialdec)
        :rtype: bool
        """
        cdef PARTIALDECOMP* cpp_partialdec = partialdec.thisptr
        return self.thisptr.addPartialdecToFinished(cpp_partialdec)

    def addPartialdecToFinishedUnchecked(self, PartialDecomposition partialdec):
        """adds a partialdec to finished partialdecs without checking for duplicates, dev has to check this on his own

        .. seealso:: * :meth:`addPartialdecToFinished()`

        :param partialdec: pointer of partialdec that is going to be added unchecked to the finished partialdecs (data structure to carry finished decompositions)
        :type partialdec: :class:`PartialDecomposition`
        """
        cdef PARTIALDECOMP* cpp_partialdec = partialdec.thisptr
        self.thisptr.addPartialdecToFinishedUnchecked(cpp_partialdec)

    def addVarPartition(self, VarPart partition):
        """adds a variable partition if it is no duplicate of an existing variable partition

        :param partition: variable partition object
        :type partition: :class:`VarPart`
        """
        cdef VarPartition* cpp_partition = partition.varPartition
        self.thisptr.addVarPartition(cpp_partition)

    def clearAncestorPartialdecs(self):
        """clears ancestor partialdec data structure

        .. note:: does not free the partialdecs themselves
        """
        self.thisptr.clearAncestorPartialdecs()

    def clearCurrentPartialdecs(self):
        """clears current partialdec data structure

        .. note:: does not free the partialdecs themselves.
        """
        self.thisptr.clearCurrentPartialdecs()

    def clearFinishedPartialdecs(self):
        """clears finished partialdec data structure

        .. note:: does not free the partialdecs themselves.
        """
        self.thisptr.clearFinishedPartialdecs()

    def createConssAdjacency(self):
        """creates the constraint adjacency datastructure that is used (if created) for some methods to faster access the constarints that have variables in common
        """
        self.thisptr.createConssAdjacency()

    def freeTemporaryData(self):
        """frees temporary data that is only needed during the detection process
        """
        self.thisptr.freeTemporaryData()

    def getAncestorPartialdec(self, partialdecindex):
        """returns a partialdec from ancestor partialdec data structure with given index

        :param partialdecindex: index of the partialdec
        :type partialdecindex: int
        :return: partialdec from ancestor partialdec
        :rtype: :class:`PartialDecomposition`
        """
        cdef PARTIALDECOMP* result = self.thisptr.getAncestorPartialdec(partialdecindex)
        return PartialDecomposition.create(result)

    def getConsPartition(self, partitionIndex):
        """returns a constraint partition

        :param partitionIndex: index of a constraint partition
        :type partitionIndex: int
        :return: constraint partition object
        :rtype: :class:`ConsPart`
        """
        cdef ConsPartition* result = self.thisptr.getConsPartition(partitionIndex)
        return ConsPart.create(result, <DetProbData>weakref.proxy(self))

    def getConsPartitions(self):
        """returns list to stored constraint partitions

        :return: list to stored constraint partitions
        :rtype: list
        """
        cdef nconsparts = self.thisptr.getNConsPartitions()
        return [ConsPart.create(self.thisptr.getConsPartition(r), self) for r in range(nconsparts)]

    def getCons(self, consIndex):
        """returns constraint related to a constraint index

        :param consIndex: index of a constraint
        :type consIndex: int
        :return: constraint related to a constraint index
        :rtype: :class:`Constraint`
        """
        return Constraint.create(self.thisptr.getCons(consIndex))

    def getConssForCons(self, consIndex):
        """returns list of constraint indices that have a common variable with the given constraint

        .. note:: constraint adjacency data structure has to initilized

        :param consIndex: index of constraint
        :type consIndex: int
        :return: list of constraint indices that have a common variable with the given constraint
        :rtype: list
        """
        return self.thisptr.getConssForCons(consIndex)

    def getConssForVar(self, varIndex):
        """returns the constraint indices of the coefficient matrix for a variable

        :param varIndex: index of variable
        :type varIndex: int
        :return: list of constraint indices that have a nonzero entry with this variable
        :rtype: list
        """
        return self.thisptr.getConssForVar(varIndex)

    def getOpenPartialdecs(self):
        """determines all partialdecs from current (open) partialdec data structure

        :return: all partialdecs in current (open) partialdec data structure
        :rtype: list
        """
        cdef vector[PARTIALDECOMP*] result = self.thisptr.getOpenPartialdecs()
        return [PartialDecomposition.create(r) for r in result]

    def getFinishedPartialdec(self, int partialdecindex):
        """returns a partialdec from finished partialdec data structure

        :param partialdecindex: index of partialdec
        :type partialdecindex: int
        :return: partialdec from finished partialdec data structure
        :rtype: :class:`PartialDecomposition`
        """
        cdef PARTIALDECOMP* result = self.thisptr.getFinishedPartialdec(partialdecindex)
        return PartialDecomposition.create(result)

    def getFinishedPartialdecs(self):
        """gets all finished partialdecs

        :return: all finished partialdecs
        :rtype: list
        """
        cdef vector[PARTIALDECOMP*] result = self.thisptr.getFinishedPartialdecs()
        return [PartialDecomposition.create(r) for r in result]

    def getIndexForCons(self, Constraint cons):
        """returns the constraint index related to a constraint

        :param cons: constraint the index is asked for
        :type cons: :class:`Constraint`
        :return: the constraint index related to constraint
        :rtype: int
        """
        return self.thisptr.getIndexForCons(cons.scip_cons)

    def getIndexForVar(self, Variable var):
        """returns the variable index related to a variable

        :param var: variable the index is asked for
        :type var: :class:`Variable`
        :return: the variable index related to variable
        :rtype: int
        """
        return self.thisptr.getIndexForVar(var.scip_var)

    def getNAncestorPartialdecs(self):
        """returns size of ancestor partialdec data structure

        :return: size of ancestor partialdec data structure
        :rtype: int
        """
        return self.thisptr.getNAncestorPartialdecs()

    def getNConsPartitions(self):
        """returns number of different constraint partitions

        :return: number of different constraint partitions
        :rtype: int
        """
        return self.thisptr.getNConsPartitions()

    def getNConss(self):
        """returns the number of variables considered in the detprobdata

        :return: number of variables considered in the detprobdata
        :rtype: int
        """
        return self.thisptr.getNConss()

    def getNConssForCons(self, consIndex):
        """returns the number of constraints for a given constraint index

        :param consIndex: index of constraint
        :type consIndex: int
        :return: the number of constraints for a given constraint index
        :rtype: int
        """
        return self.thisptr.getNConssForCons(consIndex)

    def getNConssForVar(self, varIndex):
        """returns the number of constraints for a given variable index where the variable has a nonzero entry in

        :param varIndex: index of variable
        :type varIndex: int
        :return: the number of constraints for a given variable
        :rtype: int
        """
        return self.thisptr.getNConssForVar(varIndex)

    def getNOpenPartialdecs(self):
        """returns size of current (open) partialdec data structure

        :return: size of current (open) partialdec data structure
        :rtype: int
        """
        return self.thisptr.getNOpenPartialdecs()

    def getNFinishedPartialdecs(self):
        """size of finished partialdec data structure

        :return: size of finished partialdec data structure
        :rtype: int
        """
        return self.thisptr.getNFinishedPartialdecs()

    def getNPartialdecs(self):
        """returns the number of stored partialdecs

        :return: number of stored partialdecs
        :rtype: int
        """
        return self.thisptr.getNPartialdecs()

    def getNNonzeros(self):
        """returns the number of nonzero entries in the coefficient matrix

        :return: the number of nonzero entries in the coefficient matrix
        :rtype: int
        """
        return self.thisptr.getNNonzeros()

    def getNVarPartitions(self):
        """returns number of different variable partitions

        :return: number of different variable partitions
        :rtype: int
        """
        return self.thisptr.getNVarPartitions()

    def getNVars(self):
        """return the number of variables considered in the detprobdata

        :return: the number of variables considered in the detprobdata
        :rtype: int
        """
        return self.thisptr.getNVars()

    def getNVarsForCons(self, consIndex):
        """returns the number of variables for a given constraint

        :param consIndex: index of constraint
        :type consIndex: int
        :return: the number of variables for a given constraint
        :rtype: int
        """
        return self.thisptr.getNVarsForCons(consIndex)

    def getOrigVarsFixedZero(self):
        """returns list of all original variables that are fixed to zero

        :return: list of variables
        :rtype: list
        """
        cdef vector[SCIP_VAR*] result = self.thisptr.getOrigVarsFixedZero()
        return [Variable.create(v) for v in result]

    def getRelevantConss(self):
        """returns list of all constraints that are not marked as deleted or obsolete

        :return: list of constraints
        :rtype: list
        """
        cdef vector[SCIP_CONS*] result = self.thisptr.getRelevantConss()
        return [Constraint.create(c) for c in result]

    def getRelevantVars(self):
        """returns list of all problem variables that are not fixed to 0

        :return: list of variables
        :rtype: list
        """
        cdef vector[SCIP_VAR *] result = self.thisptr.getRelevantVars()
        return [Variable.create(v) for v in result]

    def getModel(self):
        """returns the corresponding Model instance wrapping the scip data structure

        :return: the corresponding Model instance wrapping scip data structure
        :rtype: :class:`pyscipopt.Model`
        """
        cdef SCIP* scip = self.thisptr.getScip()
        return SCIPModel.create(scip)

    def getSortedCandidatesNBlocks(self, object candidates):
        """gets the candidates for number of blocks added by the user followed by the found ones sorted in descending order by how often a candidate was proposed

        :param candidates: will contain the candidates for number of blocks sorted in descending order by how often a candidate was added
        :type candidates: object
        """
        cdef vector[int] cpp_candidates = candidates
        self.thisptr.getSortedCandidatesNBlocks(cpp_candidates)

    def getVal(self, row, col):
        """returns a coefficient from the coefficient matrix

        :param row: index of the constraint to be considered
        :type row: int
        :param col: index of the variable to be considered
        :type col: int
        :return: a coefficient from the coefficient matrix
        :rtype: double
        """
        return self.thisptr.getVal(row, col)

    def getValsForCons(self, consIndex):
        """returns the nonzero coefficients of the coefficient matrix for a constraint

        .. note:: same order as in :meth:`getVarsForCons`

        :param consIndex: index of constraint
        :type consIndex: int
        :return: list of coefficients of in matrix for constraints
        :rtype: list
        """
        return self.thisptr.getValsForCons(consIndex)

    def getVarPartition(self, partitionIndex):
        """returns variable partition with given index

        :param partitionIndex: index of variable partition
        :type partitionIndex: int
        :return: variable partition with given index
        :rtype: :class:`VarPart`
        """
        cdef VarPartition* result = self.thisptr.getVarPartition(partitionIndex)
        return VarPart.create(result, self)

    def getVarPartitions(self):
        """returns list to stored variable partitions

        :return: list to stored variable partitions
        :rtype: list
        """
        cdef vector[VarPartition*] result = self.thisptr.getVarPartitions()
        return [VarPart.create(r, self) for r in result]

    def getVar(self, varIndex):
        """returns variable related to a variable index

        :return: variable related to a variable index
        :rtype: :class:`Variable`
        """
        return Variable.create(self.thisptr.getVar(varIndex))

    def getVarsForCons(self, consIndex):
        """returns the variable indices of the coefficient matrix for a constraint

        :param consIndex: index of constraint
        :type consIndex: int
        :return: list of the variable indices of the coefficient matrix for a constraint
        :rtype: list
        """
        return self.thisptr.getVarsForCons(consIndex)

    def isConsCardinalityCons(self, consIndex):
        """returns whether a constraint is a cardinality constraint, i.e. of the .. math::`\\sum_{i} x_i = b`

        :param consIndex: index of constraint
        :type consIndex: int
        :return: returns whether a constraint is a cardinality constraint
        :rtype: bool
        """
        return self.thisptr.isConsCardinalityCons(consIndex)

    def isConssAdjInitialized(self):
        """determines whether or not the constraint-constraint adjacency data structure is initilized

        :return: True iff the constraint-constraint adjacency data structure is initilized
        :rtype: bool
        """
        return self.thisptr.isConssAdjInitialized()

    def isConsSetpp(self, consIndex):
        """is cons with specified indec partitioning, or packing covering constraint?

        :param consIndex: index of constraint
        :type consIndex: int
        :return: is constraint with specified index partitioning or packing covering constraint
        :rtype: bool
        """
        return self.thisptr.isConsSetpp(consIndex)

    def isConsSetppc(self, consIndex):
        """is cons with specified index partitioning packing, or covering constraint?

        :param consIndex: index of constraint
        :type consIndex: int
        :return: whether a constraint is partitioning packing, or covering constraint?
        :rtype: bool
        """
        return self.thisptr.isConsSetppc(consIndex)

    def isPartialdecDuplicateofFinished(self, PartialDecomposition partialdec):
        """checks if partialdec is a duplicate of an existing finished partialdec

        :param partialdec: partialdec to be checked
        :type partialdec: :class:`PartialDecomposition`
        :return: True iff partialdec is a duplicate of an existing finished partialdec
        :rtype: bool
        """
        cdef PARTIALDECOMP* cpp_partialdec = partialdec.thisptr
        return self.thisptr.isPartialdecDuplicateofFinished(cpp_partialdec)

    def isAssignedToOrigProb(self):
        """returns True if the matrix structure corresponds to the presolved problem

        :return: True if the matrix structure corresponds to the presolved problem
        :rtype: bool
        """
        return self.thisptr.isAssignedToOrigProb()

    def partialdecIsNoDuplicateOfPartialdecs(self, PartialDecomposition comppartialdec, object partialdecs, sort):
        """check if partialdec is a duplicate of any given partialdecs

        :param comppartialdec: partialdec to be checked
        :type comppartialdec: :class:`PartialDecomposition`
        :param partialdecs: partialdecs to compare with
        :type partialdecs: object
        :param sort: sort the vars and conss data structures in the partialdecs by their indices
        :type sort: bool
        :return: True iff partialdec is no duplicate of any given partialdecs
        :rtype: bool
        """
        cdef PARTIALDECOMP* cpp_comppartialdec = comppartialdec.thisptr
        cdef vector[PARTIALDECOMP*] cpp_partialdecs
        cdef PARTIALDECOMP* partialdecs_ptr = NULL
        cdef PartialDecomposition partialdecs_element
        for partialdecs_element in partialdecs:
            partialdecs_ptr = <PARTIALDECOMP*> partialdecs_element.thisptr
            cpp_partialdecs.push_back(partialdecs_ptr)
        return self.thisptr.partialdecIsNoDuplicateOfPartialdecs(cpp_comppartialdec, cpp_partialdecs, sort)

    def sortFinishedForScore(self):
        """sorts partialdecs in finished partialdecs data structure according to the current scoretype
        """
        self.thisptr.sortFinishedForScore()

    def translatePartialdecs(self, DetProbData otherdata, object otherpartialdecs):
        """translates partialdecs if the index structure of the problem has changed, e.g. due to presolving

        :param otherdata: old detprobdata
        :type otherdata: :class:`DetProbData`
        :param otherpartialdecs: partialdecs to be translated
        :type otherpartialdecs: object
        :return: translated partialdecs
        :rtype: list
        """
        cdef DETPROBDATA* cpp_otherdata = otherdata.thisptr
        cdef vector[PARTIALDECOMP*] cpp_otherpartialdecs
        cdef PARTIALDECOMP* otherpartialdecs_ptr = NULL
        cdef PartialDecomposition otherpartialdecs_element
        for otherpartialdecs_element in otherpartialdecs:
            otherpartialdecs_ptr = <PARTIALDECOMP*> otherpartialdecs_element.thisptr
            cpp_otherpartialdecs.push_back(otherpartialdecs_ptr)
        cdef vector[PARTIALDECOMP*] result = self.thisptr.translatePartialdecs(cpp_otherdata, cpp_otherpartialdecs)
        return [PartialDecomposition.create(r) for r in result]
