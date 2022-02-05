cdef class DetProbData:
    """class to manage the detection process and data for one coefficient matrix of a MIP, usually there is one detprobdata for the original and one detprobdata for the presolved problem.
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
    #     @param _originalProblem True iff the detprobdata is created for the presolved problem.
    #     """
    #     # TODO implement function
    #     raise NotImplementedError()

    @property
    def candidatesNBlocks(DetProbData self):
        """candidate for the number of blocks, second int indicates how often a candidate was added.
        """
        cdef vector[pair[int, int]] result = self.thisptr.candidatesNBlocks
        return result

    @candidatesNBlocks.setter
    def candidatesNBlocks(DetProbData self, object candidatesNBlocks):
        cdef vector[pair[int, int]] cpp_candidatesNBlocks = candidatesNBlocks
        self.thisptr.candidatesNBlocks = cpp_candidatesNBlocks

    @property
    def conspartitioncollection(DetProbData self):
        """collection of different constraint class distributions.
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
        """collection of different variable class distributions.
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
        """time that was consumed by the classification of the constraint and variables classifiers.
        """
        cdef double result = self.thisptr.classificationtime
        return result

    @classificationtime.setter
    def classificationtime(DetProbData self, double classificationtime):
        cdef double cpp_classificationtime = classificationtime
        self.thisptr.classificationtime = cpp_classificationtime

    @property
    def nblockscandidatescalctime(DetProbData self):
        """time that was used to calulate the candidates of te block number.
        """
        cdef double result = self.thisptr.nblockscandidatescalctime
        return result

    @nblockscandidatescalctime.setter
    def nblockscandidatescalctime(DetProbData self, double nblockscandidatescalctime):
        cdef double cpp_nblockscandidatescalctime = nblockscandidatescalctime
        self.thisptr.nblockscandidatescalctime = cpp_nblockscandidatescalctime

    @property
    def postprocessingtime(DetProbData self):
        """time that was spent in postproceesing decomposigtions.
        """
        cdef double result = self.thisptr.postprocessingtime
        return result

    @postprocessingtime.setter
    def postprocessingtime(DetProbData self, double postprocessingtime):
        cdef double cpp_postprocessingtime = postprocessingtime
        self.thisptr.postprocessingtime = cpp_postprocessingtime

    @property
    def translatingtime(DetProbData self):
        """time that was spent by transforming partialdecs between presolved and orig problem.
        """
        cdef double result = self.thisptr.translatingtime
        return result

    @translatingtime.setter
    def translatingtime(DetProbData self, double translatingtime):
        cdef double cpp_translatingtime = translatingtime
        self.thisptr.translatingtime = cpp_translatingtime

    def addConsPartition(DetProbData self, ConsPart partition):
        """adds a constraint partition if it is no duplicate of an existing constraint partition.
        """
        cdef ConsPartition * cpp_partition = partition.thisptr
        self.thisptr.addConsPartition(cpp_partition)

    def addCandidatesNBlocksNVotes(DetProbData self, int candidate, int nvotes):
        """adds a candidate for block number and counts how often a candidate is added.
        """
        cdef int cpp_candidate = candidate
        cdef int cpp_nvotes = nvotes
        self.thisptr.addCandidatesNBlocksNVotes(cpp_candidate, cpp_nvotes)

    def addPartialdecToAncestor(DetProbData self, PartialDecomposition partialdec):
        """adds a partialdec to ancestor partialdecs

        :param partialdec: partialdec that is added to the ancestor partialdecs.
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        self.thisptr.addPartialdecToAncestor(cpp_partialdec)

    def addPartialdecToOpen(DetProbData self, PartialDecomposition partialdec):
        """adds a partialdec to current partialdecs (data structure for partialdecs that are goin to processed in the propagation rounds)

        :param partialdec: pointer of partialdec to be added
        :return: True if the partialdecs was successfully added (i.e. it is no duplicate of a known partialdec)
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        cdef bool result = self.thisptr.addPartialdecToOpen(cpp_partialdec)
        return result

    def addPartialdecToFinished(DetProbData self, PartialDecomposition partialdec):
        """adds a partialdec to finished partialdecs

        :param partialdec: pointer of partialdec that is going to be added to the finished partialdecs (data structure to carry finished decompositions)
        :return: True if the partialdecs was successfully added (i.e. it is no duplicate of a known partialdec)

        .. seealso:: * :meth:`addPartialdecToFinishedUnchecked()`
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        cdef bool result = self.thisptr.addPartialdecToFinished(cpp_partialdec)
        return result

    def addPartialdecToFinishedUnchecked(DetProbData self, PartialDecomposition partialdec):
        """adds a partialdec to finished partialdecs without checking for duplicates, dev has to check this on his own

        :param partialdec: pointer of partialdec that is going to be added unchecked to the finished partialdecs (data structure to carry finished decompositions)

        .. seealso:: * :meth:`addPartialdecToFinished()`
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        self.thisptr.addPartialdecToFinishedUnchecked(cpp_partialdec)

    def addVarPartition(DetProbData self, VarPart partition):
        """adds a variable partition if it is no duplicate of an existing variable partition

        :param partition: varpartition to be added.
        """
        cdef VarPartition * cpp_partition = partition.thisptr
        self.thisptr.addVarPartition(cpp_partition)

    def clearAncestorPartialdecs(DetProbData self):
        """clears ancestor partialdec data structure,
        .. note:: does not free the partialdecs themselves.
        """
        self.thisptr.clearAncestorPartialdecs()

    def clearCurrentPartialdecs(DetProbData self):
        """clears current partialdec data structure

        .. note:: does not free the partialdecs themselves.
        """
        self.thisptr.clearCurrentPartialdecs()

    def clearFinishedPartialdecs(DetProbData self):
        """clears finished partialdec data structure

        .. note:: does not free the partialdecs themselves.
        """
        self.thisptr.clearFinishedPartialdecs()

    def createConssAdjacency(DetProbData self):
        """create the constraint adjacency datastructure that is used (if created) for some methods to faster access the constarints that have variables in common.
        """
        self.thisptr.createConssAdjacency()

    def freeTemporaryData(DetProbData self):
        """frees temporary data that is only needed during the detection process.
        """
        self.thisptr.freeTemporaryData()

    def getAncestorPartialdec(DetProbData self, int partialdecindex):
        """returns a partialdec from ancestor partialdec data structure with given index

        :return: partialdec from ancestor partialdec data structure.
        """
        cdef int cpp_partialdecindex = partialdecindex
        cdef PARTIALDECOMP * result = self.thisptr.getAncestorPartialdec(cpp_partialdecindex)
        return PartialDecomposition.create(result)

    def getConsPartition(DetProbData self, int partitionIndex):
        """returns pointer to a constraint partition

        :return: pointer to a cosntraint partition with the given index.
        """
        cdef int cpp_partitionIndex = partitionIndex
        cdef ConsPartition * result = self.thisptr.getConsPartition(cpp_partitionIndex)
        return ConsPart.create(result, <DetProbData>weakref.proxy(self))

    def getCons(DetProbData self, int consIndex):
        """returns the SCIP constraint related to a constraint index

        :return: the SCIP constraint related to a constraint index.
        """
        return Constraint.create(self.thisptr.getCons(consIndex))

    def getConssForCons(DetProbData self, int consIndex):
        """return array of constraint indices that have a common variable with the given constraint

        :return: return vector of constraint indices that have a common variable with the given constraint

        .. note:: constraint adjacency data structure has to initilized.
        """
        cdef int cpp_consIndex = consIndex
        cdef vector[int] result = self.thisptr.getConssForCons(cpp_consIndex)
        return result

    def getConssForVar(DetProbData self, int varIndex):
        """returns the constraint indices of the coefficient matrix for a variable

        :return: vector of constraint indices that have a nonzero entry with this variable.
        """
        cdef int cpp_varIndex = varIndex
        cdef vector[int] result = self.thisptr.getConssForVar(cpp_varIndex)
        return result

    def getOpenPartialdecs(DetProbData self):
        """determines all partialdecs from current (open) partialdec data structure

        :return:  all partialdecs in current (open) partialdec data structure
        """
        cdef vector[PARTIALDECOMP *] result = self.thisptr.getOpenPartialdecs()
        return [PartialDecomposition.create(r) for r in result]

    def getFinishedPartialdec(DetProbData self, int partialdecindex):
        """returns a partialdec from finished partialdec data structure

        :return:  partialdec from finished partialdec data structure.
        """
        cdef int cpp_partialdecindex = partialdecindex
        cdef PARTIALDECOMP * result = self.thisptr.getFinishedPartialdec(cpp_partialdecindex)
        return PartialDecomposition.create(result)

    def getFinishedPartialdecs(DetProbData self):
        """gets all finished partialdecs

        :return: all finished partialdecs.
        """
        cdef vector[PARTIALDECOMP *] result = self.thisptr.getFinishedPartialdecs()
        return [PartialDecomposition.create(r) for r in result]

    def getIndexForCons(DetProbData self, Constraint cons):
        """returns the constraint index related to a SCIP constraint

        :param cons: the SCIP constraint pointer the index is asked for
        :return: the constraint index related to a SCIP constraint.
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
        """returns size of ancestor partialdec data structure

        :return: size of ancestor partialdec data structure.
        """
        cdef int result = self.thisptr.getNAncestorPartialdecs()
        return result

    def getNConsPartitions(DetProbData self):
        """returns number of different constraint partitions

        :return: number of different constraint partitions.
        """
        cdef int result = self.thisptr.getNConsPartitions()
        return result

    def getNConss(DetProbData self):
        """returns the number of variables considered in the detprobdata

        :return: number of variables considered in the detprobdata.
        """
        cdef int result = self.thisptr.getNConss()
        return result

    def getNConssForCons(DetProbData self, int consIndex):
        """returns the number of constraints for a given constraint

        :return: the number of constraints for a given constraint.
        """
        cdef int cpp_consIndex = consIndex
        cdef int result = self.thisptr.getNConssForCons(cpp_consIndex)
        return result

    def getNConssForVar(DetProbData self, int varIndex):
        """returns the number of constraints for a given variable where the var has a nonzero entry in

        :return: the number of constraints for a given variable.
        """
        cdef int cpp_varIndex = varIndex
        cdef int result = self.thisptr.getNConssForVar(cpp_varIndex)
        return result

    def getNOpenPartialdecs(DetProbData self):
        """returns size of current (open) partialdec data structure

        :return: size of current (open) partialdec data structure.
        """
        cdef int result = self.thisptr.getNOpenPartialdecs()
        return result

    def getNFinishedPartialdecs(DetProbData self):
        """size of finished partialdec data structure

        :return:  size of finished partialdec data structure.
        """
        cdef int result = self.thisptr.getNFinishedPartialdecs()
        return result

    def getNPartialdecs(DetProbData self):
        """returns the number of stored partialdecs

        :return:  number of stored partialdecs.
        """
        cdef int result = self.thisptr.getNPartialdecs()
        return result

    def getNNonzeros(DetProbData self):
        """returns the number of nonzero entries in the coefficient matrix

        :return: the number of nonzero entries in the coefficient matrix.
        """
        cdef int result = self.thisptr.getNNonzeros()
        return result

    def getNVarPartitions(DetProbData self):
        """returns number of different variable partitions

        :return:  number of different variable partitions.
        """
        cdef int result = self.thisptr.getNVarPartitions()
        return result

    def getNVars(DetProbData self):
        """return the number of variables considered in the detprobdata

        :return: the number of variables considered in the detprobdata.
        """
        cdef int result = self.thisptr.getNVars()
        return result

    def getNVarsForCons(DetProbData self, int consIndex):
        """returns the number of variables for a given constraint

        :return: the number of variables for a given constraint.
        """
        cdef int cpp_consIndex = consIndex
        cdef int result = self.thisptr.getNVarsForCons(cpp_consIndex)
        return result

    def getOrigVarsFixedZero(DetProbData self):
        """returns pointers to all orig vars that are fixed to zero

        :return: vector of vars.
        """
        cdef vector[SCIP_VAR *] result = self.thisptr.getOrigVarsFixedZero()
        return [Variable.create(v) for v in result]

    def getRelevantConss(DetProbData self):
        """returns pointers to all constraints that are not marked as deleted or obsolete

        :return: vector of conss.
        """
        cdef vector[SCIP_CONS *] result = self.thisptr.getRelevantConss()
        return [Constraint.create(c) for c in result]

    def getRelevantVars(DetProbData self):
        """returns pointers to all problem vars that are not fixed to 0

        :return: vector of vars.
        """
        cdef vector[SCIP_VAR *] result = self.thisptr.getRelevantVars()
        return [Variable.create(v) for v in result]

    def getModel(DetProbData self):
        """returns the corresponding Model instance wrapping the scip data structure

        :return: the corresponding Model instance wrapping scip data structure.
        """
        cdef SCIP * scip = self.thisptr.getScip()
        return SCIPModel.create(scip)

    def getSortedCandidatesNBlocks(DetProbData self, object candidates):
        """gets the candidates for number of blocks added by the user followed by the found ones sorted in descending order by how often a candidate was proposed

        :param candidates: will contain the candidates for number of blocks sorted in descending order by how often a candidate was added.
        """
        cdef vector[int] cpp_candidates = candidates
        self.thisptr.getSortedCandidatesNBlocks(cpp_candidates)

    def getVal(DetProbData self, int row, int col):
        """returns a coefficient from the coefficient matrix

        :return: a coefficient from the coefficient matrix.
        """
        cdef int cpp_row = row
        cdef int cpp_col = col
        cdef double result = self.thisptr.getVal(cpp_row, cpp_col)
        return result

    def getValsForCons(DetProbData self, int consIndex):
        """returns the nonzero coefficients of the coefficient matrix for a constraint

        :return: vector of coefficients of in matrix for constraints

        :note: same order as in :meth:`getVarsForCons`.
        """
        cdef int cpp_consIndex = consIndex
        cdef vector[double] result = self.thisptr.getValsForCons(cpp_consIndex)
        return result

    def getVarPartition(DetProbData self, int partitionIndex):
        """returns pointer to a variable partition with given index

        :return: pointer to a variable partition with given index.
        """
        cdef int cpp_partitionIndex = partitionIndex
        cdef VarPartition * result = self.thisptr.getVarPartition(cpp_partitionIndex)
        return VarPart.create(result)

    def getVarPartitions(DetProbData self):
        """returns vector to stored variable partitions

        :return: returns vector to stored variable partitions.
        """
        cdef vector[VarPartition *] result = self.thisptr.getVarPartitions()
        return [VarPart.create(r) for r in result]

    def getVar(DetProbData self, int varIndex):
        """returns SCIP variable related to a variable index

        :return: SCIP variable pointer related to a variable index.
        """
        cdef int cpp_varIndex = varIndex
        # TODO implement function
        raise NotImplementedError()

    def getVarsForCons(DetProbData self, int consIndex):
        """returns the variable indices of the coefficient matrix for a constraint

        :return: the variable indices of the coefficient matrix for a constraint.
        """
        cdef int cpp_consIndex = consIndex
        cdef vector[int] result = self.thisptr.getVarsForCons(cpp_consIndex)
        return result

    def isConsCardinalityCons(DetProbData self, int consindexd):
        """returns whether a constraint is a cardinality constraint, i.e. of the .. math::`\\sum_{i} x_i = b`

        :param consindexd: index of constraint that is be checked
        :return: returns whether a constraint is a cardinality constraint
        """
        cdef int cpp_consindexd = consindexd
        cdef bool result = self.thisptr.isConsCardinalityCons(cpp_consindexd)
        return result

    def isConssAdjInitialized(DetProbData self):
        """determines whether or not the constraint-constraint adjacency data structure is initilized

        :return: True iff the constraint-constraint adjacency data structure is initilized.
        """
        cdef unsigned int result = self.thisptr.isConssAdjInitialized()
        return result

    def isConsSetpp(DetProbData self, int consindexd):
        """is cons with specified indec partitioning, or packing covering constraint?

        :param consindexd: index of the given cons
        :return: is cons with specified indec partitioning, or packing covering constraint.
        """
        cdef int cpp_consindexd = consindexd
        cdef bool result = self.thisptr.isConsSetpp(cpp_consindexd)
        return result

    def isConsSetppc(DetProbData self, int consindexd):
        """is cons with specified index partitioning packing, or covering constraint?

        :param consindexd: index of cons to be checked
        :return: whether a constraint is partitioning packing, or covering constraint?.
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
        """check if partialdec is a duplicate of an existing finished partialdec

        :param partialdec: partialdec to be checked
        :return: True iff partialdec is a duplicate of an existing finished partialdec.
        """
        cdef PARTIALDECOMP * cpp_partialdec = partialdec.thisptr
        cdef unsigned int result = self.thisptr.isPartialdecDuplicateofFinished(cpp_partialdec)
        return result

    def isAssignedToOrigProb(DetProbData self):
        """returns True if the matrix structure corresponds to the presolved problem

        :return: True if the matrix structure corresponds to the presolved problem.
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
        """check if partialdec is a duplicate of any given partialdecs

        :param comppartialdec: partialdec to be checked
        :param partialdecs: partialdecs to compare with
        :param sort: sort the vars and conss data structures in the partialdecs by their indices
        :return: True iff partialdec is no duplicate of any given partialdecs.
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
        """sorts partialdecs in finished partialdecs data structure according to the current scoretype.
        """
        self.thisptr.sortFinishedForScore()

    def translatePartialdecs(DetProbData self, DetProbData otherdata, object otherpartialdecs):
        """translates partialdecs if the index structure of the problem has changed, e.g. due to presolving

        :return: translated partialdecs
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
