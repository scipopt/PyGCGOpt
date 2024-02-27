__version__ = '0.3.1'

# required for Python 3.8 on Windows
import os
if hasattr(os, 'add_dll_directory'):
    if os.getenv('SCIPOPTDIR'):
        os.add_dll_directory(os.path.join(os.getenv('SCIPOPTDIR').strip('"'), 'bin'))

# Expose pyscipopt object through pygcgopt
from pyscipopt    import multidict
from pyscipopt    import Benders
from pyscipopt    import Benderscut
from pyscipopt    import Branchrule
from pyscipopt    import Nodesel
from pyscipopt    import Conshdlr
from pyscipopt    import Eventhdlr
from pyscipopt    import Heur
from pyscipopt    import Presol
from pyscipopt    import Pricer
from pyscipopt    import Prop
from pyscipopt    import Sepa
from pyscipopt    import LP
from pyscipopt    import Expr
from pyscipopt    import quicksum
from pyscipopt    import quickprod
from pyscipopt    import exp
from pyscipopt    import log
from pyscipopt    import sqrt
from pyscipopt    import SCIP_RESULT
from pyscipopt    import SCIP_PARAMSETTING
from pyscipopt    import SCIP_PARAMEMPHASIS
from pyscipopt    import SCIP_STATUS
from pyscipopt    import SCIP_STAGE
from pyscipopt    import SCIP_PROPTIMING
from pyscipopt    import SCIP_PRESOLTIMING
from pyscipopt    import SCIP_HEURTIMING
from pyscipopt    import SCIP_EVENTTYPE
from pyscipopt    import SCIP_LPSOLSTAT
from pyscipopt    import SCIP_BRANCHDIR
from pyscipopt    import SCIP_BENDERSENFOTYPE
from pyscipopt    import SCIP_ROWORIGINTYPE

# export user-relevant objects
from pygcgopt.gcg import Model
from pygcgopt.gcg import Detector
from pygcgopt.gcg import PricingSolver
from pygcgopt.gcg import ConsClassifier
from pygcgopt.gcg import VarClassifier
from pygcgopt.gcg import Score
from pygcgopt.gcg import PY_GCG_PRICINGSTATUS as GCG_PRICINGSTATUS
from pygcgopt.gcg import PY_CONS_DECOMPINFO as CONS_DECOMPINFO
from pygcgopt.gcg import PY_VAR_DECOMPINFO as VAR_DECOMPINFO
from pygcgopt.gcg import PY_USERGIVEN as USERGIVEN
