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
        :param name: name of the class
        :type name: str
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

cdef class VarPart:

    @staticmethod
    cdef create(VarPartition* varPartition, DetProbData detProbData):
        if varPartition == NULL:
            raise Warning("cannot create VarPart with VarPartition* == NULL")
        new_VarPart = VarPart()
        new_VarPart.varPartition = varPartition
        new_VarPart.detProbData = detProbData
        return new_VarPart

    def getClassDescription(self, classindex):
        """returns the description of the class corresponding to the classindex

        :param classindex: index of the class
        :type classindex: int
        :return: description of the class
        :rtype: str
        """
        return self.varPartition.getClassDescription(classindex).decode('utf-8')

    def getClassName(self, classindex):
        """returns the name of the class corresponding to the classindex

        :param classindex: index of the class
        :type classindex: int
        :return: name of the class
        :rtype: str
        """
        return self.varPartition.getClassName(classindex).decode('utf-8')

    def getName(self):
        """returns the name of the variable partition

        :return: name of the constraint partition
        :rtype: str
        """
        return self.varPartition.getName().decode('utf-8')

    def getNClasses(self):
        """returns the number of classes of the variable partition

        :return: number of classes
        :rtype: int
        """
        return self.varPartition.getNClasses()

    def removeEmptyClasses(self):
        """removes all classes which do not have any assigned indices (classindices may change)

        :return: number of removed classes
        :rtype: int
        """
        return self.varPartition.removeEmptyClasses()

    def setClassDescription(self, classindex, desc):
        """sets the description of the class corresponding to the classindex

        :param classindex: index of the class
        :type classindex: int
        :param desc: description of the class
        :type desc: str
        """
        c_desc = str_conversion(desc)
        self.varPartition.setClassDescription(classindex, c_desc)

    def setClassName(self, classindex, name):
        """sets the name of the class corresponding to the classindex

        :param classindex: index of the class
        :type classindex: int
        :param name: name of the class
        :type name: str
        """
        c_name = str_conversion(name)
        self.varPartition.setClassName(classindex, name)

    def addClass(self, name, desc, decompInfo):
        """adds a new class to the variable partition

        :param name: name of the class
        :type name: str
        :param desc: description of the class
        :type desc: str
        :param decompInfo:                          #missing
        :type decompInfo: :class:`VAR_DECOMPINFO`
        :return: index of the added class
        :rtype: int
        """
        c_name = str_conversion(name)
        c_desc = str_conversion(desc)
        return self.varPartition.addClass(c_name, c_desc, decompInfo)

    def assignVarToClass(self, Variable var, classindex):
        """assigns a variable to the class corresponding to the classindex

        :param var: variable that should be assigned to a class
        :type var: :class:`Variable`
        :param classindex: index of the class
        :type classindex: int
        """
        cdef int cpp_varindex = self.detProbData.getIndexForVar(var)
        self.varPartition.assignVarToClass(cpp_varindex, classindex)

    def getAllSubsets(self, all, linking, master, block): #change resulting varindex in list to Variable
        """returns a vector containing all possible subsets of the chosen classindices

        :param all:
        :type all: bool
        :param linking:
        :type linking: bool
        :param master:
        :type master: bool
        :param block:
        :type block: bool
        :return: list of lists containing all 
        :rtype: list
        """
        return self.varPartition.getAllSubsets(all, linking, master, block)

    def getClassDecompInfo(self, classindex):
        return self.varPartition.getClassDecompInfo(classindex)

    def getClassNameOfVar(self, Variable var):
        cdef int cpp_varindex = self.detProbData.getIndexForVar(var)
        return self.varPartition.getClassNameOfVar(cpp_varindex).decode('utf-8')

    def getClassOfVar(self, Variable var):
        cdef int cpp_varindex = self.detProbData.getIndexForVar(var)
        return self.varPartition.getClassOfVar(cpp_varindex)

    def getNVars(self):
        """returns the number of variables of the variable partition

        :return: number of variables
        :rtype: int
        """
        return self.varPartition.getNVars()

    def getNVarsOfClasses(self):
        return self.varPartition.getNVarsOfClasses()

    def isVarClassified(self, Variable var):
        cdef int cpp_varindex = self.detProbData.getIndexForVar(var)
        return self.varPartition.isVarClassified(cpp_varindex)

    def reduceClasses(self, maxNumberOfClasses):
        cdef VarPartition* result = self.varPartition.reduceClasses(maxNumberOfClasses)
        return VarPart.create(result, self.detProbData)

    def setClassDecompInfo(self, classindex, decompInfo):
        self.varPartition.setClassDecompInfo(classindex, decompInfo)
