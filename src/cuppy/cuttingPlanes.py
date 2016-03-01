__version__    = '0.5.2'
__author__     = 'Aykut Bulut and Ted Ralphs'
__license__    = 'Eclipse Public License'
__maintainer__ = 'Ted Ralphs'
__email__      = 'ted@lehigh.edu'
__url__        = 'https://github.com/tkralphs/CuPPy'

import sys
import math
import importlib as ilib
import numpy as np
from cylp.cy import CyClpSimplex
from cylp.py.modeling import CyLPArray
PYOMO_INSTALLED = True
try:
    from pyomo.environ import *
except ImportError:
    PYOMO_INSTALLED = False

DISPLAY_ENABLED = True
try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    DISPLAY_ENABLED = False

def isInt(x, epsilon):
    '''
    Return True if x is an integer, or if x is a numpy array
    with all integer elements, False otherwise
    '''
    if isinstance(x, (int, long, float)):
        return abs(math.floor(x) - x) < epsilon
    return (np.abs(np.around(x) - x) < epsilon).all()

def getFraction(x):
    'Return the fraction part of x: x - floor(x)'
    return x - math.floor(x)
    
def gomoryCut(lp, integerIndices = None, sense = '>=', sol = None,
              rowInds = None, value = None, epsilon = .01):
    '''Return the Gomory cut of rows in ``rowInds`` of lp 
    (a CyClpSimplex object)'''
    cuts = []
    if sol == None:
        sol = lp.primalVariableSolution['x']
    if rowInds is None:
        rowInds = range(lp.nConstraints)
    if integerIndices is None:
        integerIndices = range(lp.nVariables)
    for row in rowInds:
        basicVarInd = lp.basicVariables[row]
        if (basicVarInd in integerIndices) and (not isInt(sol[basicVarInd], epsilon)):
            f0 = getFraction(sol[basicVarInd])
            f = []
            for i in range(lp.nVariables):
                if i in lp.basicVariables:
                    #This is to try to avoid getting very small numbers that 
                    #should be zero
                    f.append(0)
                else:
                    f.append(getFraction(lp.tableau[row, i]))
            pi = np.array([f[j]/f0 if f[j] <= f0 
                           else (1-f[j])/(1-f0) for j in range(lp.nVariables)])
            pi_slacks = np.array([x/f0 if x > 0 else -x/(1-f0)  
                                 for x in lp.tableau[row, lp.nVariables:]])
            pi -= pi_slacks * lp.coefMatrix
            pi0 = (1 - np.dot(pi_slacks, lp.constraintsUpper) if sense == '<='
                   else 1 + np.dot(pi_slacks, lp.constraintsUpper))
            if sense == '>=':
                cuts.append((pi, pi0))
            else:
                cuts.append((-pi, -pi0))
    return cuts
            
def disp_relaxation(A, b, cuts = [], sol = None):
    #TODO: Check sense of inequalities by looking explicitly at
    #      lp.constraintsUpper and lp.constraintsLower
    p = Polyhedron2D(A = A, b = b)
    f = Figure()
    f.add_polyhedron(p, label = 'Polyhedron $P$')
    f.set_xlim([p.xlim[0], p.xlim[1]])
    f.set_ylim([p.ylim[0], p.ylim[1]])
    pI = p.make_integer_hull()
    f.add_polyhedron(pI, show_int_points = True, color = 'red',
                     linestyle = 'dashed',
                     label = 'Convex hull of integer points')
    for (coeff, r) in cuts:
        f.add_line(coeff, r, p.xlim, p.ylim, color = 'green', linestyle = 'dashed')
    if sol is not None:
        f.add_point(sol, radius = .05)
    f.show()

def read_instance(module_name = None, file_name = None):

    if module_name is not None:
        lp = CyClpSimplex()

        mip = ilib.import_module(module_name)
            
        A = np.matrix(mip.A)
        #print np.linalg.cond(A)
        b = CyLPArray(mip.b)
        
        #Warning: At the moment, you must put bound constraints in explicitly for split cuts
        x_l = CyLPArray([0 for _ in range(mip.numVars)])
            
        x = lp.addVariable('x', mip.numVars)
        
        lp += x >= x_l
        try:
            x_u = CyLPArray(getattr(mip, 'x_u'))
            lp += x <= x_u
        except:
            x_u = None        
        lp += (A * x <= b if mip.sense[1] == '<=' else
               A * x >= b)
        c = CyLPArray(mip.c)
        lp.objective = -c * x if mip.sense[0] == 'Max' else c * x
        return lp, x, mip.A, mip.b, mip.sense[1], mip.integerIndices
    elif file_name is not None:
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
        lpNew = CyClpSimplex()
        x = lpNew.addVariable('x', lp.nCols)
        lpNew += A * x <= b
#        Trying to force bounds to be part of the A matrix
#        Doesn't work
#        for i in range(lp.nCols):
#            e = [0 for _ in range(lp.nCols)]
#            row = [0 for _ in range(lp.nCols)]
#            if lp.variablesUpper[i] < infinity:
#                e[i] = 1
#                lpNew += CyLPArray(e) * x <= lp.variablesUpper[i]
#            if lp.variablesLower[i] > -infinity:
#                e[i] = -1
#                lpNew += CyLPArray(e) * x <= -lp.variablesLower[i]
        lpNew += x <= lp.variablesUpper
        lpNew += x >= lp.variablesLower
        lpNew.objective = lp.objective
        return lpNew, x, None, None, '<=', integerIndices
    else:
        print "No file or module name specified..."
        return None, None, None, None, None, None

def solve(module_name = None, file_name = None, whichCuts = [], debug_print = False,
          epsilon = .01, max_iter = 100, max_cuts = 10, display = False):    

    if not DISPLAY_ENABLED:
        display = False
    
    if module_name is not None:
        lp, x, A, b, sense, integerIndices = read_instance(module_name = module_name)
    elif file_name is not None:
        lp, x, A, b, sense, integerIndices = read_instance(file_name = file_name)
    else:
        print "No file or module name specified...exiting"
        
    if lp.nCols > 2:
        display = False
    lp.logLevel = 0
    
    if display:
        disp_relaxation(A, b)
    
    for i in xrange(max_iter):
        print 'Iteration ', i
        lp.primal(startFinishOptions = 'x')
        print 'Current bound:', lp.objectiveValue
        #Binv = np.zeros(shape = (lp.nConstraints, lp.nConstraints))
        #for i in range(lp.nVariables, lp.nVariables+lp.nConstraints):
        #    lp.getBInvACol(i, Binv[i-lp.nVariables,:])
        #rhs = lp.rhs
        if sense == '<=':
            rhs = np.dot(lp.basisInverse, lp.constraintsUpper)
        else:
            rhs = np.dot(lp.basisInverse, lp.constraintsLower)
        sol = lp.primalVariableSolution['x']
        if debug_print:
            print 'Current basis inverse:'
            print lp.basisInverse
            print 'Condition number of basis inverse'
            print np.linalg.cond(lp.basisInverse)
            print "Current tableaux:"
            print lp.tableau
            print "Current right hand side:\n", rhs
            print lp.rhs
        print 'Current solution: ', sol
        if isInt(sol[integerIndices], epsilon):
            print 'Integer solution found!'
            break
        cuts = []
        for (cg, args) in whichCuts:
            cuts += cg(lp, integerIndices, sense, sol, **args)
        if cuts == []:
            print 'No cuts found!'
            break
        if display:
            disp_relaxation(A, b, cuts, sol)
        for (coeff, r) in cuts[:max_cuts]:
            #TODO sort cuts by degree of violation
            if sense == '<=':
                print 'Adding cut: ', coeff, '<=', r
                lp += CyLPArray(coeff) * x <= r
            else:
                print 'Adding cut: ', coeff, '>=', r
                lp += CyLPArray(coeff) * x >= r
            if display:
                A.append(coeff.tolist())
                b.append(r)
    
    if display:
        disp_relaxation(A, b)

if __name__ == '__main__':
    
    solve(module_name = 'MIP6', whichCuts = [(gomoryCut, {})], display = True)




