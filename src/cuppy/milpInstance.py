'''
Created on Apr 26, 2016

@author: TedR
'''
from builtins import range
from builtins import object
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

class MILPInstance(object):
    def __init__(self, module_name = None, file_name = None,
                 A = None, b = None, c = None,
                 points = None, rays = None,
                 l = None, u = None,
                 sense = None, integerIndices = None, 
                 numVars = None):
        
        if file_name is not None:
            # We got a file name, so ignore everything else and read in the instance
            lp = CyClpSimplex()
            lp.extractCyLPModel(file_name)
            self.integerIndices = [i for (i, j) in enumerate(lp.integerInformation) if j == True]
            infinity = lp.getCoinInfinity()
            self.A = lp.coefMatrix
            self.b = CyLPArray([0 for _ in range(lp.nRows)])
            for i in range(lp.nRows):
                if lp.constraintsLower[i] > -infinity:
                    if lp.constraintsUpper[i] < infinity:
                        raise Exception('Cannot handle ranged constraints')
                    self.b[i] = -lp.constraintsLower[i]
                    for j in range(lp.nCols):
                        self.A[i, j] = -self.A[i, j]
                elif lp.constraintsUpper[i] < infinity:
                    self.b[i] = lp.constraintsUpper[i]
                else:
                    raise Exception('Constraint with no bounds detected')
            #x = lp.addVariable('x', lp.nCols)
            #lp += A * x <= b
            #lp += x <= lp.variablesUpper
            #lp += x >= lp.variablesLower
            #lp.objective = lp.objective
            self.c = lp.objective
            self.sense = '<='
            numVars = lp.nCols
            self.lp = lp
            self.l = lp.variablesLower
            self.u = lp.variablesUpper
            self.x = lp.getVarByName('x')
        else:
            min_or_max = None
            if module_name is not None:
                # We got a module name, read the data from there
                mip = ilib.import_module(module_name)
                self.A = mip.A if hasattr(mip, 'A') else None
                self.b = mip.b if hasattr(mip, 'b') else None
                points = mip.points if hasattr(mip, 'points') else None
                rays = mip.rays if hasattr(mip, 'rays') else None
                self.c = mip.c if hasattr(mip, 'c') else None
                self.sense = mip.sense[1] if hasattr(mip, 'sense') else None
                min_or_max = mip.sense[0] if hasattr(mip, 'sense') else None
                self.integerIndices = mip.integerIndices if hasattr(mip, 'integerIndices') else None
                self.u = CyLPArray(mip.u) if hasattr(mip, 'u') else None
                self.l = CyLPArray(mip.l) if hasattr(mip, 'l') else None
                numVars = mip.numVars if hasattr(mip, 'numVars') else None
                self.x_sep = mip.x_sep if hasattr(mip, 'x_sep') else None
                if numVars is None and mip.A is not None:
                    numVars = len(mip.A)
       
                if numVars is None:
                    raise "Must specify number of variables when problem is not"   
            else:
                self.A = A
                self.b = b
                self.c = c
                self.l = l
                self.u = u
                self.points = points
                self.rays = rays
                if sense is not None:
                    self.sense = sense[1]
                    min_or_max = sense[0]
                self.integerIndices = integerIndices
                
            lp = CyClpSimplex()
            if self.A is not None:
                A = np.matrix(self.A)
                b = CyLPArray(self.b)
            elif numVars <= 2 and GRUMPY_INSTALLED:
                p = Polyhedron2D(points = points, rays = rays)
                A = np.matrix(p.hrep.A)
                b = np.matrix(p.hrep.b)
            else:
                raise "Must specify problem in inequality form with more than two variables\n"   
        
                
            x = lp.addVariable('x', numVars)

            if l is not None:
                lp += x >= l
            if u is not None:
                lp += x <= u
            lp += (A * x <= b if self.sense == '<=' else
                   A * x >= b)
            c = CyLPArray(self.c)
            if min_or_max == 'Max':
                lp.objective = -c * x
            else:
                lp.objective = c * x
            self.lp = lp
            self.x = x
