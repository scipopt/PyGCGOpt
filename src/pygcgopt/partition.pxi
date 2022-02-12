cdef class ConsPart:
    """Base class holding a pointer to corresponding ConsPartition"""

    @staticmethod
    cdef create(ConsPartition* consPartition, DetProbData detProbData):
        if consPartition == NULL:
            raise Warning("cannot create ConsPart with ConsPartition* == NULL")
        new_ConsPart = ConsPart()
        new_ConsPart.consPartition = consPartition
        new_ConsPart.detProbData = detProbData
        return new_ConsPart

    def getClassDescription(self, classindex):
        """returns the description of the class corresponding to the classindex

        :param classindex: index of the class
        :type classindex: int
        :return: description of the class
        :rtype: str
        """
        return self.consPartition.getClassDescription(classindex).decode('utf-8')

    def getClassName(self, classindex):
        """returns the name of the class corresponding to the classindex

        :param classindex: index of the class
        :type classindex: int
        :return: name of the class
        :rtype: str
        """
        return self.consPartition.getClassName(classindex).decode('utf-8')

    def getName(self):
        """returns the name of the constraint partition

        :return: name of the constraint partition
        :rtype: str
        """
        return self.consPartition.getName().decode('utf-8')

    def getNClasses(self):
        """returns the number of classes of the constraint partition

        :return: number of classes
        :rtype: int
        """
        return self.consPartition.getNClasses()

    def removeEmptyClasses(self):
        """removes all classes which do not have any assigned indices (classindices may change)

        :return: number of removed classes
        :rtype: int
        """
        return self.consPartition.removeEmptyClasses()

    def setClassDescription(self, classindex, desc):
        """sets the description of the class corresponding to the classindex

        :param classindex: index of the class
        :type classindex: int
        :param desc: description of the class
        :type desc: str
        """
        c_desc = str_conversion(desc)
        self.consPartition.setClassDescription(classindex, c_desc)

    def setClassName(self, classindex, name):
        """sets the name of the class corresponding to the classindex

        :param classindex: index of the class
        :type classindex: int
        :param desc: name of the class
        :type desc: str
        """
        c_name = str_conversion(name)
        self.consPartition.setClassName(classindex, name)

    def addClass(self, name, desc, decompInfo):
        c_name = str_conversion(name)
        c_desc = str_conversion(desc)
        return self.consPartition.addClass(c_name, c_desc, decompInfo)

    def assignConsToClass(self, Constraint cons, classindex):
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        self.consPartition.assignConsToClass(cpp_consindex, classindex)

    def getAllSubsets(self, both, only_master, only_pricing):
        return self.consPartition.getAllSubsets(both, only_master, only_pricing)

    def getClassDecompInfo(self, classindex):
        return self.consPartition.getClassDecompInfo(classindex)

    def getClassNameOfCons(self, Constraint cons):
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        return self.consPartition.getClassNameOfCons(cpp_consindex).decode('utf-8')

    def getClassOfCons(self, Constraint cons):
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        return self.consPartition.getClassOfCons(cpp_consindex)

    def getNConss(self):
        return self.consPartition.getNConss()

    def getNConssOfClasses(self):
        return self.consPartition.getNConssOfClasses()

    def isConsClassified(self, Constraint cons):
        cdef int cpp_consindex = self.detProbData.getIndexForCons(cons)
        return self.consPartition.isConsClassified(cpp_consindex)

    def reduceClasses(self, maxNumberOfClasses):
        cdef ConsPartition* result = self.consPartition.reduceClasses(maxNumberOfClasses)
        return ConsPart.create(result, self.detProbData)

    def setClassDecompInfo(self, classindex, decompInfo):
        self.consPartition.setClassDecompInfo(classindex, decompInfo)

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

    def setClassName(VarPart self, classindex, name):
        """sets the name of the class
        """
        c_name = str_conversion(name)
        self.thisptr.setClassName(classindex, c_name)

    def setClassDescription(VarPart self, classindex, desc):
        """sets the description of the class
        """
        c_desc = str_conversion(desc)
        self.thisptr.setClassDescription(classindex, c_desc)

    def setClassDecompInfo(VarPart self, classindex, decompInfo):
        """sets the decomposition code of a class

        :param classindex: index of the class
        :type classindex: int
        :param decompInfo: decomposition code of class
        :type decompInfo: VAR_DECOMPINFO
        """
        self.thisptr.setClassDecompInfo(classindex, decompInfo)

    def removeEmptyClasses(VarPart self):
        """removes all classes which do not have any assigned index (classindices may change)
        """
        return self.thisptr.removeEmptyClasses()

    def getName(VarPart self):
        """returns the name of the partition
        """
        return self.thisptr.getName().decode('utf-8')

    def getNClasses(VarPart self):
        """returns the number of classes the partition provides
        """
        return self.thisptr.getNClasses()
