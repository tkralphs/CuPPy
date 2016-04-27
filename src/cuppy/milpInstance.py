'''
Created on Apr 26, 2016

@author: TedR
'''
from cylp.cy import CyClpSimplex
from cylp.py.modeling import CyLPArray
import importlib as ilib
import numpy as np
GRUMPY_INSTALLED = True
try:
    from src.grumpy.polyhedron2D import Polyhedron2D
except ImportError:
    try:
        from coinor.grumpy.polyhedron2D import Polyhedron2D
    except ImportError:
        GRUMPY_INSTALLED = False

class MILPInstance:
    def __init__(self, module_name = None, file_name = None,
                 A = None, b = None, c = None,
                 points = None, rays = None,
                 sense = None, integerIndices = None, 
                 numVars = None):

        if file_name is not None:
            # We got a file name, so ignore everything else and read in the instance
            lp = CyClpSimplex()
            lp.extractCyLPModel(file_name)
            integerIndices = [i for (i, j) in enumerate(lp.integerInformation) if j == True]
            infinity = lp.getCoinInfinity()
            A = lp.coefMatrix
            b = CyLPArray([0 for _ in range(lp.nRows)])
            for i in range(lp.nRows):
                if lp.constraintsLower[i] > -infinity:
                    if lp.constraintsUpper[i] < infinity:
                        raise Exception('Cannot handle ranged constraints')
                    b[i] = -lp.constraintsLower[i]
                    for j in range(lp.nCols):
                        A[i, j] = -A[i, j]
                elif lp.constraintsUpper[i] < infinity:
                    b[i] = lp.constraintsUpper[i]
                else:
                    raise Exception('Constraint with no bounds detected')
            x = lp.addVariable('x', lp.nCols)
            lp += A * x <= b
            lp += x <= lp.variablesUpper
            lp += x >= lp.variablesLower
            lp.objective = lp.objective
    
        else:
            if module_name is not None:
                # We got a module name, read the data from there
                mip = ilib.import_module(module_name)
                A = mip.A if hasattr(mip, 'A') else None
                b = mip.b if hasattr(mip, 'b') else None
                c = mip.c if hasattr(mip, 'c') else None
                points = mip.points if hasattr(mip, 'points') else None
                rays = mip.rays if hasattr(mip, 'rays') else None
                sense = mip.sense if hasattr(mip, 'sense') else None
                numVars = mip.numVars if hasattr(mip, 'numVars') else None
                self.x_sep = mip.x_sep if hasattr(mip, 'x_sep') else None
       
            lp = CyClpSimplex()
            if A is not None:
                A = np.matrix(A)
                b = CyLPArray(b)
            elif numVars <= 2 and GRUMPY_INSTALLED:
                p = Polyhedron2D(points = points, rays = rays)
                A = np.matrix(p.hrep.A)
                b = np.matrix(p.hrep.b)
            else:
                raise("Must specify problem in inequality form with more than two variables\n"+
                      "or when GrUMPy is not installed")   
        
            #Warning: At the moment, you must put bound constraints in explicitly for split cuts
            x_l = CyLPArray([0 for _ in range(numVars)])
                
            x = lp.addVariable('x', numVars)
            
            lp += x >= x_l
            try:
                x_u = CyLPArray(getattr(mip, 'x_u'))
                lp += x <= x_u
            except:
                x_u = None        
            lp += (A * x <= b if sense[1] == '<=' else
                   A * x >= b)
            c = CyLPArray(c)
            lp.objective = -c * x if sense[0] == 'Max' else c * x
            self.lp = lp
            self.x = x
            self.A = A
            self.b = b
            self.sense = sense[1]
            self.integerIndices = integerIndices
