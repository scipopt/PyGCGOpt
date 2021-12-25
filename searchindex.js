Search.setIndex({docnames:["CHANGELOG","INSTALL","README","api/decomposition","api/detector","api/detprobdata","api/gcgmodel","api/index","api/pricingsolver","examples/alldecomps/alldecomps","examples/cpmp/cpmp","examples/index","index"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":4,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":3,"sphinx.domains.rst":2,"sphinx.domains.std":2,nbsphinx:3,sphinx:56},filenames:["CHANGELOG.md","INSTALL.md","README.md","api/decomposition.rst","api/detector.rst","api/detprobdata.rst","api/gcgmodel.rst","api/index.rst","api/pricingsolver.rst","examples/alldecomps/alldecomps.ipynb","examples/cpmp/cpmp.ipynb","examples/index.rst","index.md"],objects:{"pygcgopt.gcg":[[5,0,1,"","DetProbData"],[4,0,1,"","Detector"],[6,0,1,"","GCGModel"],[3,0,1,"","PartialDecomposition"],[8,0,1,"","PricingSolver"]],"pygcgopt.gcg.DetProbData":[[5,1,1,"","addCandidatesNBlocksNVotes"],[5,1,1,"","addConsPartition"],[5,1,1,"","addPartialdecToAncestor"],[5,1,1,"","addPartialdecToFinished"],[5,1,1,"","addPartialdecToFinishedUnchecked"],[5,1,1,"","addPartialdecToOpen"],[5,1,1,"","addVarPartition"],[5,2,1,"","candidatesNBlocks"],[5,2,1,"","classificationtime"],[5,1,1,"","clearAncestorPartialdecs"],[5,1,1,"","clearCurrentPartialdecs"],[5,1,1,"","clearFinishedPartialdecs"],[5,2,1,"","conspartitioncollection"],[5,1,1,"","createConssAdjacency"],[5,1,1,"","freeTemporaryData"],[5,1,1,"","getAncestorPartialdec"],[5,1,1,"","getCons"],[5,1,1,"","getConsPartition"],[5,1,1,"","getConssForCons"],[5,1,1,"","getConssForVar"],[5,1,1,"","getFinishedPartialdec"],[5,1,1,"","getFinishedPartialdecs"],[5,1,1,"","getIndexForCons"],[5,1,1,"","getModel"],[5,1,1,"","getNAncestorPartialdecs"],[5,1,1,"","getNConsPartitions"],[5,1,1,"","getNConss"],[5,1,1,"","getNConssForCons"],[5,1,1,"","getNConssForVar"],[5,1,1,"","getNFinishedPartialdecs"],[5,1,1,"","getNNonzeros"],[5,1,1,"","getNOpenPartialdecs"],[5,1,1,"","getNPartialdecs"],[5,1,1,"","getNVarPartitions"],[5,1,1,"","getNVars"],[5,1,1,"","getNVarsForCons"],[5,1,1,"","getOpenPartialdecs"],[5,1,1,"","getOrigVarsFixedZero"],[5,1,1,"","getRelevantConss"],[5,1,1,"","getRelevantVars"],[5,1,1,"","getSortedCandidatesNBlocks"],[5,1,1,"","getVal"],[5,1,1,"","getValsForCons"],[5,1,1,"","getVar"],[5,1,1,"","getVarPartition"],[5,1,1,"","getVarPartitions"],[5,1,1,"","getVarsForCons"],[5,1,1,"","isAssignedToOrigProb"],[5,1,1,"","isConsCardinalityCons"],[5,1,1,"","isConsSetpp"],[5,1,1,"","isConsSetppc"],[5,1,1,"","isConssAdjInitialized"],[5,1,1,"","isPartialdecDuplicateofFinished"],[5,2,1,"","nblockscandidatescalctime"],[5,1,1,"","partialdecIsNoDuplicateOfPartialdecs"],[5,2,1,"","postprocessingtime"],[5,1,1,"","sortFinishedForScore"],[5,1,1,"","translatePartialdecs"],[5,2,1,"","translatingtime"],[5,2,1,"","varpartitioncollection"]],"pygcgopt.gcg.Detector":[[4,2,1,"","detectorname"],[4,1,1,"","exitDetector"],[4,1,1,"","finishPartialdec"],[4,1,1,"","freeDetector"],[4,1,1,"","initDetector"],[4,2,1,"","model"],[4,1,1,"","postprocessPartialdec"],[4,1,1,"","propagatePartialdec"],[4,1,1,"","setParamAggressive"],[4,1,1,"","setParamDefault"],[4,1,1,"","setParamFast"]],"pygcgopt.gcg.GCGModel":[[6,1,1,"","addDecomposition"],[6,1,1,"","addDecompositionFromConss"],[6,1,1,"","addPreexistingPartialDecomposition"],[6,1,1,"","addVar"],[6,1,1,"","createDecomposition"],[6,1,1,"","createPartialDecomposition"],[6,1,1,"","detect"],[6,1,1,"","getDualbound"],[6,1,1,"","getMasterProb"],[6,1,1,"","includeDefaultPlugins"],[6,1,1,"","includeDetector"],[6,1,1,"","includePricingSolver"],[6,1,1,"","listDecompositions"],[6,1,1,"","listDetectors"],[6,1,1,"","listPricingSolvers"],[6,1,1,"","optimize"],[6,1,1,"","presolve"],[6,1,1,"","printStatistics"],[6,1,1,"","printVersion"],[6,1,1,"","setDetectorEnabled"],[6,1,1,"","setDetectorFinishingEnabled"],[6,1,1,"","setDetectorPostprocessingEnabled"],[6,1,1,"","setGCGSeparating"],[6,1,1,"","setPricingSolverEnabled"],[6,1,1,"","setPricingSolverExactEnabled"],[6,1,1,"","setPricingSolverHeuristicEnabled"],[6,1,1,"","writeAllDecomps"]],"pygcgopt.gcg.PartialDecomposition":[[3,1,1,"","addAncestorID"],[3,1,1,"","addBlock"],[3,1,1,"","addClockTime"],[3,1,1,"","addDecChangesFromAncestor"],[3,1,1,"","addDetectorChainInfo"],[3,1,1,"","addNNewBlocks"],[3,1,1,"","addPctConssFromFree"],[3,1,1,"","addPctConssToBlock"],[3,1,1,"","addPctConssToBorder"],[3,1,1,"","addPctVarsFromFree"],[3,1,1,"","addPctVarsToBlock"],[3,1,1,"","addPctVarsToBorder"],[3,1,1,"","aggInfoCalculated"],[3,1,1,"","alreadyAssignedConssToBlocks"],[3,1,1,"","assignCurrentStairlinking"],[3,1,1,"","assignOpenConssToMaster"],[3,1,1,"","assignPartialdecFromConstoblockVector"],[3,1,1,"","assignSmallestComponentsButOneConssAdjacency"],[3,2,1,"","bendersScore"],[3,2,1,"","borderAreaScore"],[3,1,1,"","buildDecChainString"],[3,1,1,"","calcAggregationInformation"],[3,1,1,"","calcStairlinkingVars"],[3,1,1,"","checkAllConssAssigned"],[3,1,1,"","checkConsistency"],[3,2,1,"","classicScore"],[3,1,1,"","complete"],[3,1,1,"","completeByConnected"],[3,1,1,"","completeByConnectedConssAdjacency"],[3,1,1,"","completeGreedily"],[3,1,1,"","considerImplicits"],[3,1,1,"","copy"],[3,1,1,"","copyPartitionStatistics"],[3,1,1,"","deleteEmptyBlocks"],[3,1,1,"","deleteOpencons"],[3,1,1,"","deleteOpenvar"],[3,1,1,"","displayInfo"],[3,1,1,"","findVarsLinkingToMaster"],[3,1,1,"","findVarsLinkingToStairlinking"],[3,1,1,"","fixConsToBlock"],[3,1,1,"","fixConsToBlockByName"],[3,1,1,"","fixConsToBlockId"],[3,1,1,"","fixConsToMaster"],[3,1,1,"","fixConsToMasterByName"],[3,1,1,"","fixConssToBlock"],[3,1,1,"","fixConssToBlockId"],[3,1,1,"","fixConssToMaster"],[3,1,1,"","fixVarToBlock"],[3,1,1,"","fixVarToBlockByName"],[3,1,1,"","fixVarToLinking"],[3,1,1,"","fixVarToLinkingByName"],[3,1,1,"","fixVarToMaster"],[3,1,1,"","fixVarToMasterByName"],[3,1,1,"","fixVarToStairlinking"],[3,1,1,"","getAncestorID"],[3,1,1,"","getAncestorList"],[3,1,1,"","getBlockConss"],[3,1,1,"","getBlocksForRep"],[3,1,1,"","getConssForBlocks"],[3,1,1,"","getDetectorClockTime"],[3,1,1,"","getDetectorClockTimes"],[3,1,1,"","getDetectorchainInfo"],[3,1,1,"","getDetprobdata"],[3,1,1,"","getFinishedByFinisher"],[3,1,1,"","getHashValue"],[3,1,1,"","getID"],[3,1,1,"","getLinkingvars"],[3,1,1,"","getMasterconss"],[3,1,1,"","getMastervars"],[3,1,1,"","getNAncestors"],[3,1,1,"","getNBlocks"],[3,1,1,"","getNCoeffsForBlock"],[3,1,1,"","getNCoeffsForMaster"],[3,1,1,"","getNConss"],[3,1,1,"","getNConssForBlock"],[3,1,1,"","getNDetectors"],[3,1,1,"","getNLinkingvars"],[3,1,1,"","getNMasterconss"],[3,1,1,"","getNMastervars"],[3,1,1,"","getNNewBlocks"],[3,1,1,"","getNNewBlocksVector"],[3,1,1,"","getNOpenconss"],[3,1,1,"","getNOpenvars"],[3,1,1,"","getNReps"],[3,1,1,"","getNStairlinkingvars"],[3,1,1,"","getNTotalStairlinkingvars"],[3,1,1,"","getNVars"],[3,1,1,"","getNVarsForBlock"],[3,1,1,"","getNVarsForBlocks"],[3,1,1,"","getOpenconss"],[3,1,1,"","getOpenconssVec"],[3,1,1,"","getOpenvars"],[3,1,1,"","getOpenvarsVec"],[3,1,1,"","getPctConssFromFree"],[3,1,1,"","getPctConssFromFreeVector"],[3,1,1,"","getPctConssToBlock"],[3,1,1,"","getPctConssToBlockVector"],[3,1,1,"","getPctConssToBorder"],[3,1,1,"","getPctConssToBorderVector"],[3,1,1,"","getPctVarsFromFree"],[3,1,1,"","getPctVarsFromFreeVector"],[3,1,1,"","getPctVarsToBlock"],[3,1,1,"","getPctVarsToBlockVector"],[3,1,1,"","getPctVarsToBorder"],[3,1,1,"","getPctVarsToBorderVector"],[3,1,1,"","getRepForBlock"],[3,1,1,"","getRepVarmap"],[3,1,1,"","getStairlinkingvars"],[3,1,1,"","getTranslatedpartialdecid"],[3,1,1,"","getUsergiven"],[3,1,1,"","getVarProbindexForBlock"],[3,1,1,"","getVarsForBlock"],[3,1,1,"","hasSetppMaster"],[3,1,1,"","hasSetppcMaster"],[3,1,1,"","hasSetppccardMaster"],[3,1,1,"","isAssignedToOrigProb"],[3,1,1,"","isComplete"],[3,1,1,"","isConsMastercons"],[3,1,1,"","isConsOpencons"],[3,2,1,"","isSelected"],[3,1,1,"","isTrivial"],[3,1,1,"","isVarBlockvarOfBlock"],[3,1,1,"","isVarLinkingvar"],[3,1,1,"","isVarMastervar"],[3,1,1,"","isVarOpenvar"],[3,1,1,"","isVarStairlinkingvar"],[3,1,1,"","isVarStairlinkingvarOfBlock"],[3,2,1,"","maxForWhiteAggScore"],[3,2,1,"","maxForWhiteScore"],[3,2,1,"","maxWhiteScore"],[3,2,1,"","partForWhiteAggScore"],[3,2,1,"","partForWhiteScore"],[3,1,1,"","prepare"],[3,1,1,"","refineToBlocks"],[3,1,1,"","refineToMaster"],[3,1,1,"","removeAncestorID"],[3,1,1,"","removeMastercons"],[3,1,1,"","setAncestorList"],[3,1,1,"","setConsPartitionStatistics"],[3,1,1,"","setConsToBlock"],[3,1,1,"","setConsToMaster"],[3,1,1,"","setDetectorClockTimes"],[3,1,1,"","setFinishedByFinisher"],[3,1,1,"","setFinishedByFinisherOrig"],[3,1,1,"","setNBlocks"],[3,1,1,"","setPctConssFromFreeVector"],[3,1,1,"","setPctConssToBlockVector"],[3,1,1,"","setPctConssToBorderVector"],[3,1,1,"","setPctVarsFromFreeVector"],[3,1,1,"","setPctVarsToBlockVector"],[3,1,1,"","setPctVarsToBorderVector"],[3,1,1,"","setStemsFromOrig"],[3,1,1,"","setTranslatedpartialdecid"],[3,1,1,"","setUsergiven"],[3,1,1,"","setVarPartitionStatistics"],[3,1,1,"","setVarToBlock"],[3,1,1,"","setVarToLinking"],[3,1,1,"","setVarToMaster"],[3,1,1,"","setVarToStairlinking"],[3,1,1,"","shouldCompletedByConsToMaster"],[3,1,1,"","showVisualization"],[3,1,1,"","sort"],[3,2,1,"","strongDecompScore"]],"pygcgopt.gcg.PricingSolver":[[8,1,1,"","exitSolution"],[8,1,1,"","exitSolver"],[8,1,1,"","freeSolver"],[8,1,1,"","initSolution"],[8,1,1,"","initSolver"],[8,2,1,"","model"],[8,1,1,"","solve"],[8,1,1,"","solveHeuristic"],[8,2,1,"","solvername"],[8,1,1,"","updateSolver"]]},objnames:{"0":["py","class","Python class"],"1":["py","method","Python method"],"2":["py","attribute","Python attribute"]},objtypes:{"0":"py:class","1":"py:method","2":"py:attribute"},terms:{"0":[1,3,5,6,9,10],"00":[9,10],"000":[9,10],"0000":10,"000000":9,"00000000000000e":9,"000000e":[9,10],"001":[9,10],"002":9,"003":9,"005":10,"009706853038245034":10,"01":[9,10],"019291161956034086":10,"02":[9,10],"024809e":10,"03":[9,10],"04":[9,10],"04f":10,"05":[9,10],"050000e":10,"051111e":10,"055556":9,"06":[9,10],"065000e":10,"07":9,"083333":9,"09":[9,10],"0s":[9,10],"1":[1,3,6,9,11],"10":[9,10],"100":[9,10],"101":10,"10167":10,"102":10,"102000e":10,"102222e":10,"10229":10,"1030k":9,"1032k":9,"10437":10,"11":[9,10],"110000e":10,"111111":9,"1154k":9,"1155k":9,"1156k":9,"1157k":9,"1158k":9,"1159k":9,"116":10,"12":[9,10],"1214k":9,"123333e":10,"1257k":9,"1258k":9,"1259k":9,"1260k":9,"1261k":9,"1262k":9,"1263k":9,"1264k":9,"1265k":9,"1266k":9,"1267k":9,"1268k":9,"1269k":9,"1270k":9,"1271k":9,"1272k":9,"1273k":9,"1274k":9,"1275k":9,"1276k":9,"1277k":9,"1278k":9,"1279k":9,"1281k":9,"1282k":9,"1284k":9,"1286k":9,"1287k":9,"13":10,"13000000000000e":10,"130000e":10,"138889":9,"14":[9,10],"140":9,"140000e":10,"1426k":9,"1427k":9,"15":[9,10],"150":10,"150000e":10,"166667":9,"1699":10,"1749":10,"1799":10,"18":[9,10],"1839":10,"1af84b0617":10,"1s":[9,10],"2":[1,3,9,11],"20":10,"200":9,"2002":10,"2010":10,"2021":10,"2053":10,"2081":10,"21":10,"22":[9,10],"222222":9,"2237":10,"23":10,"2302":10,"2320":10,"2450":10,"245000e":10,"25":[9,10],"2500":10,"2507":10,"2550":10,"28":9,"280":10,"28m":10,"29":10,"2933":10,"2983":10,"2em":10,"2s":[9,10],"3":[1,9,11],"30":10,"303":10,"3033":10,"305":10,"3055":10,"3083":10,"31":10,"3169":10,"327":10,"33":[9,10],"33333333333333e":9,"333333e":9,"3386":10,"33m":10,"34":[9,10],"340000e":10,"347":10,"34m":10,"35":10,"356":10,"357":10,"36":10,"366":10,"37":[9,10],"371":10,"380000e":10,"38m":10,"39":10,"3939":10,"395":10,"3s":[9,10],"4":[1,9,10],"40":[9,10],"41":9,"41m":10,"42":10,"43":10,"4384":10,"4482":10,"46":[9,10],"4851":10,"48514851485148514":10,"485149":10,"4853426519122501":10,"49":10,"49504950495049505":10,"4s":10,"5":[1,9,10],"50":[9,10],"500000e":9,"51":10,"5187":10,"51m":10,"52":10,"52m":10,"54m":10,"55":9,"56":[9,10],"57m":10,"580000e":10,"5815":10,"59m":10,"5s":10,"6":[1,9,10],"60":[9,10],"60m":10,"61m":10,"62m":10,"66":9,"66666666666667e":9,"666667e":9,"6679":10,"67":9,"6em":10,"6s":10,"7":[1,9,10],"711000e":10,"73":[9,10],"75":9,"76":10,"785000e":10,"7s":10,"8":[1,9,10],"848000e":10,"89":10,"8s":10,"9":[1,9,10],"928k":9,"929k":9,"95":10,"9636":10,"97":10,"9s":10,"boolean":6,"break":9,"byte":10,"case":[1,3,9,11],"class":[3,4,5,6,8,9,10],"default":[1,4,6],"do":[3,10],"export":1,"final":3,"function":[2,3,9,10,12],"import":[9,10],"int":[3,5,9,10],"new":[3,6,9,10],"public":0,"return":[3,5,6,9,10],"short":[1,3],"static":[3,9],"super":9,"true":[3,5,6,9,10],"try":10,"var":[3,5,9,10],"while":9,A:[3,6,9,10],As:[1,10],At:10,Be:1,For:[1,3,6,9,10],If:[1,6,10],In:[1,6,9,10],It:9,On:1,The:[1,3,6,11],Then:[1,9,10],There:[9,10],To:[1,3,9,10],With:[9,10],_:[9,10],__init__:9,__pycache__:1,_idx_:9,aachen:10,abort:[6,9,10],abortpricingint:9,about:3,abov:[1,3,9,10],access:[3,5,6,10],accord:[3,5,9,10],accordingli:3,accumul:9,activ:1,ad:[2,3,5,6,9,10,12],add:[1,3,5,6,9,10],addancestorid:3,addblock:3,addcandidatesnblocksnvot:5,addclocktim:3,addcon:10,addconspartit:5,adddecchangesfromancestor:3,adddecomposit:6,adddecompositionfromconss:6,adddetectorchaininfo:3,addit:[1,6,9,10],additionalnblock:3,addnnewblock:3,addpartialdectoancestor:5,addpartialdectofinish:5,addpartialdectofinisheduncheck:5,addpartialdectoopen:5,addpctconssfromfre:3,addpctconsstoblock:3,addpctconsstobord:3,addpctvarsfromfre:3,addpctvarstoblock:3,addpctvarstobord:3,addpreexistingpartialdecomposit:[6,10],address:3,addvar:[6,10],addvarpartit:5,adjac:[3,5],ado:9,after:[1,4,6,9,10],afterward:[6,9],agginfocalcul:3,aggreg:[3,10],aggregat:3,aggress:4,aka:3,al:9,alignat:10,all:[1,3,5,6,9,10],all_conss:9,all_decomps_detector:9,alldecomposit:6,alldecomps_inst:9,alldecompsdetector:9,allow:[3,10],almost:10,along:[9,10],alreadi:[3,9],alreadyassignedconsstoblock:3,also:[2,10,12],altern:3,alwai:[6,9,10],an:[1,2,3,5,6,9,10,12],ancestor:[3,5],ancestorid:3,ancestorindex:3,ani:[3,5,6,9,10],anoth:10,anymor:3,api:3,appear:3,append:[9,10],appli:[3,9,10],approach:10,appropri:[3,6,9],apt:1,ar:[1,3,5,6,9,10],area:3,arg:6,argument:[6,9],arrai:[3,5],as_posix:9,ascend:9,ask:[3,5],assert:9,assgin:3,assign:[3,9,10],assigncurrentstairlink:3,assignhittingopenconss:3,assignhittingopenvar:3,assignopenconsstomast:3,assignopenpartialhittingtomast:3,assignpartialdecfromconstoblock:3,assignpartialdecfromconstoblockvector:3,assignsmallestcomponentsbutoneconssadjac:3,automat:[3,6,11],avail:[1,10],avoid:9,awai:6,awar:1,b54569ac6:10,b:[5,10],back:3,base:[4,6,8,9,10],basi:9,becaus:9,becom:3,been:9,befor:[1,3,6,9,10],begin:10,below:9,bender:3,bendersscor:3,berlin:10,best:6,between:[5,10],bin:[1,9,10],bliss:3,block1:3,block2:3,block:[3,5,6,9,10],block_conss:6,block_id:3,blockcon:3,blockconss:3,blockid:3,blockrepid:3,blockrepres:3,blockvar:3,boilerpl:9,bookkeep:3,bool:[3,5,6],border:3,borderareascor:[3,10],bound:[6,10,11],branch:10,brand:1,brief:3,buffer:3,build:11,build_model:10,builddecchainstr:3,built:[3,9,10],c:[9,10],calcaggregationinform:3,calcstairlinkingvar:3,calcul:[3,9,10],call:[3,4,6,9,10],callback:9,calul:5,can:[1,2,3,6,9,10,12],candid:5,candidatesnblock:5,cannot:[1,3],capac:10,capacit:11,cardin:[3,5],carri:[3,5],ccon:3,cd:1,cell:10,certain:3,chain:[3,9],chang:[3,5,9,10],changelog:[2,12],check:[1,3,5],checkallconssassign:3,checkconsist:3,chg:[9,10],chosen:[4,9,10],classic:[3,10],classicscor:[3,10],classif:[5,9,10],classifi:[3,5,10],classificationtim:5,clear:5,clearancestorpartialdec:5,clearcurrentpartialdec:5,clearfinishedpartialdec:5,cliqu:[9,10],clock:3,clocktim:3,close:3,clq:[9,10],cmake:1,code:[3,9,10],coeff:[3,9,10],coeffici:[5,9,10],col:[3,5],collect:[5,9],color:9,coloring3:9,combin:9,come:9,command:[1,9,10],common:[3,5],compact:10,compar:5,compat:[1,2,12],compil:6,compl:9,complet:[1,3,4,6,10],completebyconnect:3,completebyconnectedconssadjac:3,completed_constomast:3,completegreedili:3,compon:3,comppartialdec:5,comput:[3,9],con:[3,5,9],conf:9,connect:3,connected:3,connectedbas:10,cons_pmedian:10,consadjac:3,consclass:10,consclassesmast:3,consclassifi:[9,10],consecut:3,consid:[3,5,9,10],considerimplicit:3,consindex:5,consindexd:5,consist:[3,10],consnam:3,consol:6,conspart:[3,5],conspartit:[3,9,10],conspartitioncollect:5,conss:[3,5,9,10],conss_assign:10,conss_capac:10,conss_mast:10,conss_powerset:9,conss_reform:10,consschang:8,constarint:5,constoblock:3,constomast:3,constraint:[3,5,6,9,10],constructor:9,constyp:[9,10],consult:9,consum:5,cont:[9,10],contain:[1,3,5,9],conveni:[3,6],copi:[3,9],copypartitionstatist:3,copyright:[6,10],correspond:[1,3,5],cosntraint:5,cost:[9,10],count:5,cover:[3,5,9],cpmp:10,creat:[1,3,5,6,9,10],createconssadjac:5,createdecomposit:6,createdirectori:6,createpartialdecomposit:[6,10],current:[1,3,5,6,10],current_timestamp:9,custom:11,cut:[6,9,10],cycl:9,d:10,d_:10,dantzig:[9,10],data:[3,5,9],datastructur:5,datetim:9,dcmake_build_typ:1,dcmake_install_prefix:1,dec:[6,10],dec_consclass:10,dec_densemasterconss:10,dec_neighborhoodmast:10,decchar:6,decid:[3,6,10],decincludedetector:6,decinfo:3,decomp:[9,10],decompos:[3,11],decomposigt:5,decomposit:[4,5,6,7,9,11],deeper:10,def:[9,10],defin:9,deg:[9,10],del:[9,10],delet:[3,5,9,10],deleteemptyblock:3,deleteopencon:3,deleteopenvar:3,demand:10,depend:1,deriv:9,desc:6,descend:5,destructor:4,det:9,detail:[1,3,9],detaillevel:3,detect:[3,4,5,6,9,10],detector:[3,6,7,10,11],detector_nam:6,detectorchain:3,detectorchainindex:3,detectornam:[4,6],determin:5,detprobdata:[3,4,7,9],dev:[1,5],develop:1,dgcg_dev_build:1,diabl:6,diagon:9,differ:[1,3,5,9,11],dig:10,dipopt:1,directli:[3,6,9,10],directori:[1,6],dirti:10,disabl:[6,9,10],disk:6,displai:[3,10],displayinfo:3,distanc:10,distribut:[5,10],document:[9,10],doe:[3,5,6],don:9,done:[3,9],dot:6,doubl:3,down:9,download:[1,9,10],dpapilo:1,dshare:1,dt:9,dual:[6,10,11],dual_bound:9,dualbound:[9,10],dualsolconv:8,due:[5,6,9],dump:9,duplic:5,dure:[5,6,9,10],dw:11,dyld_fallback_library_path:1,dzimpl:1,e:[1,3,5,10],each:[3,6,9,10],easili:10,either:[1,3],element:9,emphasi:4,empti:[3,6,9],enabl:[6,9],enabledfinish:6,enabledpostprocess:6,end:[9,10],enough:[6,9],ensur:[3,6],entri:5,enumer:10,environ:1,equal:6,equival:9,error:1,et:9,evalu:[3,9],everi:[3,9,10],exact:[6,10],exacten:6,exactli:[3,10],exampl:[1,9,10],exce:10,execut:[6,9,10],exhaust:[9,10],exist:[3,5,6],exist_ok:9,exit:[4,9],exitdetector:4,exitsolut:8,exitsolv:8,experi:11,explan:6,explicit:9,exploit:10,explor:[3,11],extend:10,extens:6,f:[9,10],fact:10,fals:[6,9],fast:[4,9,10],faster:5,featur:[9,10],few:9,field:10,file:[1,6,9,10],find:[1,9],findvarslinkingtomast:3,findvarslinkingtostairlink:3,finish:[3,5,6,9,10],finishingen:6,finishpartialdec:[3,4],finnish:[6,10],first:[1,3,9,10],firstblock:3,fix:[2,3,5,6,9,10,12],fixblockconss:6,fixconsstoblock:[3,6,9],fixconsstoblockid:[3,6],fixconsstomast:[3,6,9,10],fixconstoblock:[3,6,10],fixconstoblockbynam:3,fixconstoblockid:[3,6],fixconstomast:[3,6],fixconstomasterbynam:3,fixmasterconss:6,fixvartoblock:3,fixvartoblockbynam:3,fixvartolink:3,fixvartolinkingbynam:3,fixvartomast:3,fixvartomasterbynam:3,fixvartostairlink:3,fo:10,focu:3,folder:[6,10],follo:1,follow:[1,3,5,9,10],foral:10,forese:3,format:[6,10],formul:10,found:[1,2,3,5,6,9,10,12],fraction:3,free:[4,5],freedetector:4,freeprob:9,freesolv:8,freetemporarydata:5,freqcallround:6,freqcallroundorigin:6,fresh:9,from:[2,3,5,6,9,10,12],from_iter:9,fromorig:3,fuer:10,full:10,fulli:10,fun:9,further:9,furthermor:1,g:[1,5,10],gamsdomain:[9,10],gamssymbol:[9,10],gap:[9,10],gcg:[1,2,3,4,5,6,8,9,10,12],gcgconshdlrdecompcalcbendersscor:3,gcgconshdlrdecompcalcborderareascor:3,gcgconshdlrdecompcalccandidatesnblock:10,gcgconshdlrdecompcalcclassicscor:3,gcgconshdlrdecompcalcmaxforeseeingwhiteaggscor:3,gcgconshdlrdecompcalcmaxforseeingwhitescor:3,gcgconshdlrdecompcalcmaxwhitescor:3,gcgconshdlrdecompcalcsetpartforwhiteaggscor:3,gcgconshdlrdecompcalcstrongdecompositionscor:3,gcgmodel:[7,9,10],gcgoptdir:1,gener:[3,9],get:[1,3,5,9,10],get_simple_inst:10,getancestorid:3,getancestorlist:3,getancestorpartialdec:5,getblockconss:3,getblocksforrep:3,getcon:5,getconspartit:5,getconssforblock:3,getconssforcon:5,getconssforvar:5,getdetectorchaininfo:3,getdetectorclocktim:3,getdetprobdata:3,getdualbound:[6,9],getfinishedbyfinish:3,getfinishedpartialdec:5,gethashvalu:3,getid:3,getindexforcon:5,getlinkingvar:3,getmasterconss:3,getmasterprob:[6,9],getmastervar:3,getmodel:5,getnancestor:3,getnancestorpartialdec:5,getnblock:3,getncoeffsforblock:3,getncoeffsformast:3,getnconspartit:5,getnconss:[3,5,9],getnconssforblock:3,getnconssforcon:5,getnconssforvar:5,getndetector:3,getnfinishedpartialdec:5,getnlinkingvar:3,getnmasterconss:3,getnmastervar:3,getnnewblock:3,getnnewblocksvector:3,getnnonzero:5,getnopenconss:3,getnopenpartialdec:5,getnopenvar:3,getnpartialdec:5,getnrep:3,getnstairlinkingvar:3,getntotalstairlinkingvar:3,getnvar:[3,5],getnvarpartit:5,getnvarsforblock:3,getnvarsforcon:5,getopenconss:[3,9],getopenconssvec:3,getopenpartialdec:5,getopenvar:3,getopenvarsvec:3,getorigvarsfixedzero:5,getpctconssfromfre:3,getpctconssfromfreevector:3,getpctconsstoblock:3,getpctconsstoblockvector:3,getpctconsstobord:3,getpctconsstobordervector:3,getpctvarsfromfre:3,getpctvarsfromfreevector:3,getpctvarstoblock:3,getpctvarstoblockvector:3,getpctvarstobord:3,getpctvarstobordervector:3,getpresolvingtim:9,getreadingtim:9,getrelevantconss:5,getrelevantvar:5,getrepforblock:3,getrepvarmap:3,getsolvingtim:9,getsortedcandidatesnblock:5,getstairlinkingvar:3,getstatu:9,gettotaltim:9,gettranslatedpartialdecid:3,getusergiven:3,getval:5,getvalsforcon:5,getvar:5,getvarpartit:5,getvarprobindexforblock:3,getvarsforblock:3,getvarsforcon:5,ggcgconshdlrdecompcalcsetpartforseeingwhitescor:3,githash:10,given:[2,3,5,9,10,12],global:1,go:5,goal:9,goin:5,gp:3,graphic:10,greedili:3,group:10,gt:[9,10],gz:9,h:[1,9],ha:[1,3,5,9,10],happen:3,harmlessli:1,hash:3,hashabl:3,hashvalu:3,hassetppccardmast:3,hassetppcmast:3,hassetppmast:3,have:[1,3,5,6,9,10],header:1,help:10,henc:3,here:[1,9,10],heur:9,heuren:6,heurist:6,heuristicen:6,hi:5,hit:3,hold:[3,10],hole:[9,10],how:[3,5,9,10],howev:10,hspace:10,i:[3,5,9,10],id:3,ideal:1,ident:[3,10],identifi:3,iff:[3,5],ignor:3,ignoredetectionlimit:3,ij:10,impl:[9,10],implement:9,implic:[9,10],implicit:3,includ:[1,6,9,10],includedefaultplugin:6,includedetector:[6,9],includepricingsolv:6,inconsist:3,incorpor:3,increas:3,indec:5,index:[3,5,9],indic:[3,5],individu:9,inf:[9,10],info:3,inform:[2,3,6,9,12],informationstechnik:10,init_model:9,initdetector:4,initi:[0,4,6],initil:5,initsolut:8,initsolv:8,ins:6,insid:4,inspect:11,instanc:[5,6,9,11],instance_nam:9,instanci:[10,11],instead:10,instruct:[1,2,12],int_max:6,integ:3,integr:[9,10],intent:9,interact:[6,9,10],interfac:[1,2,9,10,12],intern:3,interrupt:9,introduc:3,investig:10,involv:3,ip:10,is_en:6,isassignedtoorigprob:[3,5],iscomplet:3,isconscardinalitycon:5,isconsmastercon:3,isconsopencon:3,isconssadjiniti:5,isconssetpp:5,isconssetppc:5,ispartialdecduplicateoffinish:5,isselect:[3,10],istrivi:3,isvarblockvarofblock:3,isvarlinkingvar:3,isvarmastervar:3,isvaropenvar:3,isvarstairlinkingvar:3,isvarstairlinkingvarofblock:3,item:10,iter:[3,6,9,10],iteration_idx:9,itertool:9,its:[1,3,6,9,10],itself:3,j4:1,j:10,joinpath:9,json:[9,10],jsonl:9,jupyt:[9,10],just:9,k:10,kei:9,kk:10,knapsack:10,know:[3,10],known:[5,10],konrad:10,kwarg:6,lambda:9,largest:3,last_decomp:9,later:[6,9,10],latest:[1,2,12],ld_library_path:1,least:[1,3],leav:3,left:[9,10],len:[9,10],leq:10,let:[9,10],lib:1,libgcg:1,librari:1,libscip:1,like:[1,10],limit:[3,9],line:[1,9],linear:[6,9,10],link:[3,9,10],linkingvar:10,linux:[1,3],list:[3,6,9,10],listdecomposit:[6,9,10],listdetector:[6,9],listpricingsolv:6,literatur:10,littl:9,load:10,local:1,locat:[1,10],log:[6,9],log_:9,log_filenam:9,log_path:9,logs_dir:9,longer:3,look:[9,10],loop:[4,6,9],lp:[6,9,10],lpi:1,lt:[9,10],m:[1,9,10],maco:1,magic:10,mai:9,main:[6,9],major:[1,2,12],make:[1,3,9],makefil:1,manag:[3,5],mani:3,manner:[9,10],manual:[1,10],map:[3,9],mark:5,master:[3,6,9,10],master_conss:[6,9],masterconss:10,masterpric:9,mastervar:[3,10],match:9,math:5,mathbb:10,matrix:[3,5,9,10],max:10,maxcallround:6,maxcallroundorigin:6,maxforwhiteaggscor:3,maxforwhitescor:[3,10],maximum:[3,10],maxwhit:[9,10],maxwhitescor:[3,10],mcon:[9,10],mcut:[9,10],md:[2,12],mdpt:[9,10],mean:3,measur:[9,10],median:11,medium:[9,10],mem:[9,10],memori:[4,10],menu:3,messag:1,method:[3,4,5,6,9,10],might:[3,9],min:10,mincallround:6,mincallroundorigin:6,minim:10,minut:9,mip:[5,10],miplib:[9,10],miss:[9,10],mkdir:[1,9],mlp:[9,10],mode:[6,11],model:[1,4,5,6,8,11],more:[3,9,10],most:[3,9,10],mp:[9,10],multipl:10,must:[1,10],mvar:[9,10],n:[9,10],n_cluster:10,n_conss:9,n_locat:10,name:[3,6,9],navig:1,nblock:[3,10],nblockscandidatescalctim:5,neccessari:[3,6,9],need:[1,3,5,9],neighborhoodmast:10,neither:3,new_dec:9,newlist:3,newpartialdec:9,newvector:3,next:[1,9],nlinkingvar:10,nlpi:1,nmasterconss:10,nmastervar:10,nnewblock:3,node:[6,9,10],nodebook:10,none:10,nonempti:3,nonzero:[3,5,9,10],note:[1,2,3,5,9,10,12],notebook:[9,10],noth:3,notic:[9,10],now:[9,10],nozeror:3,number:[3,5,6,9,10],nvote:5,obj:10,object:[3,6,9,10],obsolet:5,obtain:[6,9],obviou:3,occasion:1,occur:3,ocon:[9,10],off:[1,9],offer:[9,10],offet:9,often:5,onc:[3,9],one:[3,5,9,10],ones:5,onli:[1,2,3,5,9,10,12],onlin:10,open:[3,5,9,10],opencon:3,openvar:3,oper:[1,10],optim:[1,2,6,9,10,12],option:[1,6,10],order:[3,5,10],orig:[3,5],origin:[3,5,6,9,10],other:[3,9,10],otherdata:5,otherpartialdec:[3,5],otherwis:[3,6],our:[9,10],out:9,output:[3,6],ovar:[9,10],over:9,overal:3,overriden:9,overview:[3,10],own:[5,6],p550:10,p:[9,11],pack:[3,5],page:1,panic:9,paper:9,paramet:[3,5,6,9],parent:9,part:[9,10],partforwhiteaggscor:3,partforwhitescor:3,partial:[3,4,6,9],partialdec:[3,5,6,9,10],partialdecindex:5,partialdecisnoduplicateofpartialdec:5,partialdecomp:3,partialdecomposit:[3,5,6,9,10],particular:10,partit:[3,5],partitionindex:5,pass:[1,3,6,9],patch:3,path:[1,6,9,10],path_to_install_dir:1,pathlib:9,pct:3,pd:10,per:[3,9,10],percentag:3,perform:3,pip:1,plane:9,pleas:[1,2,9,10,12],plug:[6,9],plugin:[4,8,9],pointer:5,possibl:[3,9,10],postprocees:5,postprocess:[4,6,9,10],postprocessingen:6,postprocessingtim:5,postprocesspartialdec:4,potenti:3,power:9,powerset:9,precis:10,precondit:3,predefin:9,preexist:10,prepar:[3,11],presenc:3,present:[3,9],presolv:[3,5,6,9,10],presolving_tim:9,prevent:9,previou:10,previous:9,price:[6,8,9,10],pricing_solver_nam:6,pricingprob:8,pricingsolv:[6,7],primal:[9,10],primalbound:[9,10],print:[6,9,10],printstatist:6,printvers:[6,10],prior:3,prioriti:6,probe:[9,10],problem:[3,5,6,9,11],problem_nam:9,problem_path:9,probnr:8,procedur:[1,9],process:[5,9,10],produc:9,program:6,project:[2,12],propag:[3,5],propagatepartialdec:[4,9,10],propos:5,provid:[2,6,12],py:1,pygcgopt:[3,4,5,6,8,9,10],pypi:[2,12],pyscipopt:[4,8,9,10],pytest:1,python3:1,python:[1,2,3,9,10,12],q:10,q_i:10,q_j:10,quad:10,quicksum:10,r:[1,9,10],rais:3,rang:[9,10],reach:9,read:[1,9,11],read_instance_json:10,reading_tim:9,readprob:10,readproblem:9,reassign:3,receiv:9,recommend:1,recreat:10,redirectoutput:10,reduc:[9,10],referanc:9,referenc:9,refin:[3,4],refinetoblock:3,refinetomast:3,reform_conss:9,reform_conss_idx:9,reformul:[9,10],reformulation_constraint:9,regist:[3,9],relat:[3,5],relax:[6,9],releas:[0,1,2,12],relev:3,reli:3,remain:9,remov:[1,2,3,12],removeancestorid:3,removemastercon:3,renam:3,reorder:3,repid:3,replic:9,report:6,repres:[3,6],represen:3,request:3,requir:[2,9,12],rerun:[9,10],research:10,resolv:9,resourc:9,respect:3,restrict:6,result:[3,9],result_:9,results_dir:9,results_fil:9,retriev:6,right:6,rmp:9,root:[1,6],round:[5,9,10],row:[3,5,9],run:[1,6,10,11],runtim:1,rwth:10,s:[3,9,10],sake:9,same:[3,5,6,9,10],scale:[9,10],scip:[1,2,3,4,5,6,8,9,10,12],scip_error:3,scip_okai:3,scip_paramset:[6,9],scip_stag:9,scipoptdir:1,score:[3,9,10],scoretyp:5,sec:[9,10],second:[5,9,10],section:1,see:[1,2,3,6,9,12],seealso:3,seem:3,seen:10,select:[1,3,9,11],self:[3,4,5,6,8,9],sepa:9,separ:6,set:[1,3,4,6,11],setancestorlist:3,setboolparam:9,setconspartitionstatist:3,setconstoblock:3,setconstomast:3,setdetectorclocktim:3,setdetectoren:[6,9],setdetectorfinishingen:[6,9],setdetectorpostprocessingen:[6,9],setfinishedbyfinish:3,setfinishedbyfinisherorig:3,setgcgsepar:[6,9],setlogfil:9,setlongintparam:9,setminim:10,setnblock:3,setparamaggress:4,setparamdefault:4,setparamfast:4,setpartit:3,setpctconssfromfreevector:3,setpctconsstoblockvector:3,setpctconsstobordervector:3,setpctvarsfromfreevector:3,setpctvarstoblockvector:3,setpctvarstobordervector:3,setppc:10,setpricingsolveren:6,setpricingsolverexacten:6,setpricingsolverheuristicen:6,setstemsfromorig:3,settranslatedpartialdecid:3,setup:9,setusergiven:3,setvarpartitionstatist:3,setvartoblock:3,setvartolink:3,setvartomast:3,setvartostairlink:3,sever:9,share:1,should:[3,6],shouldcompletedbyconstomast:3,showvisu:3,side:[9,10],similar:9,simpl:9,simpli:[1,10],sinc:[9,10],size:[3,5],skip:[1,6],skipe:6,slp:[9,10],small:9,smaller:3,so:[1,3,10],solut:[6,9,10],solv:[6,8,9,10],solveheurist:8,solver:[1,2,6,8,9,10,12],solvernam:[6,8],solving_tim:9,some:[3,5],soplex:10,sort:[3,5,9],sortfinishedforscor:5,sourc:[2,12],specif:[1,10],specifi:[1,3,5,6,9,10],spent:5,stage:[9,10],staircas:3,stairlink:3,stairlinkingvar:3,standard:10,star:10,start:[9,10],statist:[3,6,9],statu:[1,3,9,10],stdout:6,stem:[3,9],step:9,still:[3,10],store:[3,5,6,9,10],strategi:3,strbr:9,strftime:9,string:[3,6],strong:3,strongdecompscor:3,structur:[3,5,6,9,10],studi:[9,10],sub:3,subclass:[6,9,10],subdirectori:1,success:[3,10],successfulli:5,sudo:1,suit:[1,2,12],sum:10,sum_:[5,10],summar:1,summari:11,support:10,suppos:3,sure:1,sy:9,symmetri:9,system:1,t:[9,10],tabl:1,take:[9,10],task:10,te:5,technic:10,temporari:5,term:10,terrifi:9,test:10,test_nam:1,text:10,th:3,than:1,thei:[1,3,10],them:[1,3,10],themselv:5,therefor:[1,9,10],thi:[1,2,3,5,6,9,10,12],those:9,three:[9,10],through:[6,9,10],thu:9,tighten:[9,10],time:[3,5,9,10],took:[9,10],total:[3,9,10],total_tim:9,transform:[3,5,6,9,10],translat:[3,5,9,10],translatepartialdec:5,translatingtim:5,trivial:[3,9,10],turn:6,tutori:[9,10],two:[3,9,10],txt:1,type:[3,9,10],typeerror:3,typic:10,u:6,uncheck:5,uncondit:9,unicod:[4,8],uniqu:3,univers:10,unix:1,unknown:9,unspecifi:3,up:11,updatesolv:8,upgd:[9,10],us:[1,2,3,5,6,9,11,12],useful:9,usefulrecal:6,useless:10,user:[5,6],usergiven:3,userpartialdec:3,usual:[1,2,3,5,12],utcnow:9,util:9,v:10,valu:[3,6,9,10],varbndschang:8,varclass:10,varclasseslink:3,varclassesmast:3,varclassifi:[9,10],varctor:3,variabl:[1,3,5,6,9,10],variableopen:3,varid:3,varindex:5,varmap:3,varnam:3,varobjschang:8,varobjv:[9,10],varobjvalsign:[9,10],varpart:[3,5],varpartit:[3,5,9,10],varpartitioncollect:5,vartoblock:3,vartolink:3,vartomast:3,vartostairlink:3,vartyp:[9,10],vector:[3,5],venv:1,version:[1,2,3,6,10,12],virtual:1,visual:[3,11],vtype:10,vv:10,wa:[3,5,9],walk:9,want:[1,9,10],we:[1,9,10],weather:6,well:[9,10],were:3,what:3,when:[3,4,6,9,10],where:[5,6,10],whether:[3,5],which:[1,2,3,9,10,12],white:[3,10],whose:3,why:9,wise:3,without:[3,5,6,9],witt:9,wolf:[9,10],wonder:9,work:[1,9,10],workonpartialdec:[4,9],would:[9,10],wrap:5,write:[6,9],writealldecomp:6,x:[1,10],x_:10,x_i:5,y:[9,10],y_:10,y_j:10,yet:[3,9,10],yield:[9,10],you:[1,3,9,10],your:[1,10],z:10,zentrum:10,zero:5,zib:10,zip:9,zuse:10},titles:["CHANGELOG","Installation","PyGCGOpt","Decomposition","Detector","DetProbData","GCGModel","API Reference","PricingSolver","Exploring DW Dual Bounds","Decomposing the Capacitated p-Median Problem","Examples","PyGCGOpt"],titleterms:{"1":[0,10],"2":10,"3":[0,10],"case":10,"new":1,The:[9,10],ad:0,api:7,automat:10,binari:1,bound:9,build:[1,10],capacit:10,chang:0,changelog:0,compil:1,custom:10,debug:1,decompos:10,decomposit:[3,10],detector:[4,9],detprobdata:5,differ:10,distribut:1,document:[2,12],dual:9,dw:9,exampl:11,experi:9,explor:[9,10],fix:0,from:1,gcgmodel:6,inform:1,inspect:10,instal:[1,2,12],instanc:10,instanci:9,median:10,mode:10,model:[9,10],p:10,prepar:9,pricingsolv:8,problem:10,pygcgopt:[1,2,12],pypi:1,pyscipopt:1,read:10,refer:7,remov:0,requir:1,run:9,scipoptsuit:1,select:10,set:[9,10],sourc:1,summari:[9,10],test:1,unreleas:0,up:[9,10],us:10,v0:0,visual:10}})