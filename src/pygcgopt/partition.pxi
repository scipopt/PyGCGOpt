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
        """assigns a constraint to a class.
        """
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        cdef int cpp_classindex = classindex
        self.thisptr.assignConsToClass(cpp_consindex, cpp_classindex)

    def getAllSubsets(ConsPart self, bool both, bool only_master, bool only_pricing):
        """returns a vector containing all possible subsets of the chosen classindices.
        """
        cdef bool cpp_both = both
        cdef bool cpp_only_master = only_master
        cdef bool cpp_only_pricing = only_pricing
        cdef vector[vector[int]] result = self.thisptr.getAllSubsets(cpp_both, cpp_only_master, cpp_only_pricing)
        return result

    # def getClassDecompInfo(ConsPart self, int classindex):
    #     """returns the decomposition info of a class.
    #     """
    #     cdef int cpp_classindex = classindex
    #     # TODO implement function
    #     raise NotImplementedError()

    def getClassNameOfCons(ConsPart self, Constraint cons):
        """returns the name of the class a constraint is assigned to.
        """
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        cdef const char * result = self.thisptr.getClassNameOfCons(cpp_consindex)
        return result.decode('utf-8')

    def getClassOfCons(ConsPart self, Constraint cons):
        """returns the index of the class a constraint is assigned to.
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
        """returns the number of constraints.
        """
        cdef int result = self.thisptr.getNConss()
        return result

    def getNConssOfClasses(ConsPart self):
        """returns a vector with the numbers of constraints that are assigned to the classes.
        """
        cdef vector[int] result = self.thisptr.getNConssOfClasses()
        return result

    def isConsClassified(ConsPart self, Constraint cons):
        """returns whether a constraint is already assigned to a class.
        """
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        cdef bool result = self.thisptr.isConsClassified(cpp_consindex)
        return result

    def reduceClasses(ConsPart self, int maxNumberOfClasses):
        """returns partition with reduced number of classes
        if the current number of classes is greater than an upper bound
        and lower than 2*(upper bound) (returns NULL otherwise).
        """
        cdef int cpp_maxNumberOfClasses = maxNumberOfClasses
        cdef ConsPartition * result = self.thisptr.reduceClasses(cpp_maxNumberOfClasses)
        return ConsPart.create(result, self.detProbData)

    def getName(ConsPart self):
        """returns the name of the partition"""
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
        """assigns a variable to a class.
        """
        cdef int cpp_varindex = varindex
        cdef int cpp_classindex = classindex
        self.thisptr.assignVarToClass(cpp_varindex, cpp_classindex)

    def getAllSubsets(VarPart self, bool all, bool linking, bool master, bool block):
        """returns a vector containing all possible subsets of the chosen classindices.
        """
        cdef bool cpp_all = all
        cdef bool cpp_linking = linking
        cdef bool cpp_master = master
        cdef bool cpp_block = block
        cdef vector[vector[int]] result = self.thisptr.getAllSubsets(cpp_all, cpp_linking, cpp_master, cpp_block)
        return result

    # def getClassDecompInfo(VarPart self, int classindex):
    #     """returns the decomposition info of a class.
    #     """
    #     cdef int cpp_classindex = classindex
    #     # TODO implement function
    #     raise NotImplementedError()

    def getClassNameOfVar(VarPart self, int varindex):
        """returns the name of the class a variable is assigned to.
        """
        cdef int cpp_varindex = varindex
        cdef const char * result = self.thisptr.getClassNameOfVar(cpp_varindex)
        return result.decode('utf-8')

    def getClassOfVar(VarPart self, int varindex):
        """returns the index of the class a variable is assigned to.
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
        """returns the number of variables.
        """
        cdef int result = self.thisptr.getNVars()
        return result

    def getNVarsOfClasses(VarPart self):
        """returns a vector with the numbers of variables that are assigned to the classes.
        """
        cdef vector[int] result = self.thisptr.getNVarsOfClasses()
        return result

    def isVarClassified(VarPart self, int varindex):
        """returns whether a variable is already assigned to a class.
        """
        cdef int cpp_varindex = varindex
        cdef bool result = self.thisptr.isVarClassified(cpp_varindex)
        return result

    def reduceClasses(VarPart self, int maxNumberOfClasses):
        """returns partition with reduced number of classes
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
