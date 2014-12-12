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
from cylp.py.modeling import CyLPArray, CyLPModel
from copy import deepcopy
sys.path.append('instances')

try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure

DISPLAY_ENABLED = True
try:
    import matplotlib.pyplot as plt
except ImportError:
    DISPLAY_ENABLED = False

def isInt(x):
    '''
    Return True if x is an integer, or if x is a numpy array
    with all integer elements, False otherwise
    '''
    if isinstance(x, (int, long, float)):
        return abs(math.floor(x) - x) < epsilon
    return (np.abs(np.floor(x) - x) < epsilon).all()

def getFraction(x):
    'Return the fraction part of x: x - floor(x)'
    return x - math.floor(x)

def splitCuts(lp, integerIndices = None, sense = '>=', sol = None,
              max_coeff = 1):
    A = lp.coefMatrix
    b = CyLPArray(lp.constraintsUpper)
    if integerIndices is None:
        integerIndices = range(lp.nVariables)
    if sol is None:
        sol = lp.primalVariableSolution['x']
    s = A*sol - b
    best = lp.getCoinInfinity()

    for theta in [0.1, 0.2, 0.3, 0.4, 0.5]:
        sp = CyLPModel()
        u = sp.addVariable('u', lp.nConstraints, isInt = False)
        v = sp.addVariable('v', lp.nConstraints, isInt = False)
        pi = sp.addVariable('pi', lp.nVariables, isInt = True)
        pi0 = sp.addVariable('pi0', 1, isInt = True)
    
        sp += pi + A.transpose()*u - A.transpose()*v == 0
        sp += pi0 + b*u - b*v == theta - 1
        if sense == '<=':
            sp += u >= 0
            sp += v >= 0
        else:
            #TODO this direction is not debugged
            # Is this all we need?
            sp += u <= 0
            sp += v <= 0
            
        sp.objective = (theta-1)*s*u - theta*s*v
        for i in xrange(lp.nVariables):
            if i in integerIndices:
                sp += -max_coeff <= pi[i] <= max_coeff
            else:
                sp[i] += pi[i] == 0

        cbcModel = CyClpSimplex(sp).getCbcModel()
        cbcModel.logLevel = -1
        #cbcModel.maximumSeconds = 5
        cbcModel.solve()
        if debug_print:
            print theta, cbcModel.objectiveValue
            print cbcModel.primalVariableSolution['pi'], cbcModel.primalVariableSolution['pi0']
        if cbcModel.objectiveValue < best:
            best = cbcModel.objectiveValue
            multu = cbcModel.primalVariableSolution['u']
            disjunction = cbcModel.primalVariableSolution['pi']
            rhs = cbcModel.primalVariableSolution['pi0']
            best_theta = theta
            
    alpha = A.transpose()*multu + best_theta*disjunction
    if sense == '<=':
        beta = np.dot(lp.constraintsUpper, multu) + best_theta*rhs
    else:
        beta = np.dot(lp.constraintsLower, multu) + best_theta*rhs
    if (abs(alpha) > 1e-6).any():
        return [(alpha, beta)] 
    else:
        return []
    
def gomoryCut(lp, integerIndices = None, sense = '>=', rowInds = None, value = None):
    'Return the Gomory cut of rows in ``rowInds`` of lp (a CyClpSimplex object)'
    cuts = []
    sol = lp.primalVariableSolution['x']
    if rowInds is None:
        rowInds = range(lp.nConstraints)
    if integerIndices is None:
        integerIndices = range(lp.nVariables)
    for row in rowInds:
        basicVarInd = lp.basicVariables[row]
        if (basicVarInd in integerIndices) and (not isInt(sol[basicVarInd])):
            f0 = getFraction(sol[basicVarInd])
            f = []
            for i in range(lp.nVariables):
                if i in lp.basicVariables:
                    #This is to try to avoid getting very small numbers that should be zero
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
    f.set_xlim(p.plot_min[0], p.plot_max[0])
    f.set_ylim(p.plot_min[1], p.plot_max[1])
    pI = p.make_integer_hull()
    f.add_polyhedron(pI, show_int_points = True, color = 'red',
                     linestyle = 'dashed',
                     label = 'Convex hull of integer points')
    for (coeff, r) in cuts:
        f.add_line(coeff, r, plot_max = p.plot_max, plot_min = p.plot_min,
                   color = 'green', linestyle = 'dashed')
    if sol is not None:
        f.add_point(sol, radius = .05)
    f.show()

def read_instance(module_name = True, file_name = None):

    if module_name:
        lp = CyClpSimplex()

        mip = ilib.import_module(module_name)
            
        A = np.matrix(mip.A)
        #print np.linalg.cond(A)
        b = CyLPArray(mip.b)
        
        #We assume variables have zero lower bounds
        x_l = CyLPArray([0 for _ in range(mip.numVars)])
            
        x = lp.addVariable('x', mip.numVars)
        
        lp += x >= x_l
        try:
            x_u = CyLPArray(getattr(mip, 'x_u'))
            lp += x <= x_u
        except:
            pass
        
        lp += (A * x <= b if mip.sense[1] == '<=' else
               A * x >= b)
        c = CyLPArray(mip.c)
        lp.objective = -c * x if mip.sense[0] == 'Max' else c * x
        
        return lp, x, mip.A, mip.b, mip.sense, mip.integerIndices
    else:
        #TODO Change sense of inequalities so they are all the same
        #     by explicitly checking lp.constraintsUpper and lp.constraintsLower
        #Warning: Reading MP not well tested 
        lp.extractCyLPModel(file_name)
        x = lp.cyLPModel.getVarByName('x')
        sense = ('Min', '>=')
        return lp, x, None, None, sense, integerIndices

if __name__ == '__main__':
    
    generate_splits = True
    generate_GMI = True
    debug_print = False
    epsilon = 0.01
    maxiter = 10
    max_cuts = 10
    display = True
    if not DISPLAY_ENABLED:
        display = False
    
    lp, x, A, b, sense, integerIndices = read_instance(module_name = 'MIP6')
    infinity = lp.getCoinInfinity()
    
    if display:
        disp_relaxation(A, b)
    
    for i in xrange(maxiter):
        print '______________________________________________________'
        print 'iteration ', i
        lp.primal(startFinishOptions = 'x')
        bv = lp.basicVariables
        #Binv = np.zeros(shape = (lp.nConstraints, lp.nConstraints))
        #for i in range(lp.nVariables, lp.nVariables+lp.nConstraints):
        #    lp.getBInvACol(i, Binv[i-lp.nVariables,:])
        #rhs = lp.rhs
        rhs = np.dot(lp.basisInverse, b)
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
        if isInt(sol[integerIndices]):
            print 'Integer solution found!'
            break
        cuts = []
        if generate_splits:
            cuts += splitCuts(lp, integerIndices, sense[1])
        if generate_GMI:
            cuts += gomoryCut(lp, integerIndices, sense[1])
        if cuts == []:
            print 'No cuts found!'
            break
        if display:
            disp_relaxation(A, b, cuts, sol)
        for (coeff, r) in cuts[:max_cuts]:
            #TODO sort cuts by degree of violation
            if sense[1] == '<=':
                print 'Adding cut: ', coeff, '<=', r
                lp += CyLPArray(coeff) * x <= r
            else:
                print 'Adding cut: ', coeff, '>=', r
                lp += CyLPArray(coeff) * x >= r
            A.append(coeff.tolist())
            b.append(r)
    
    disp_relaxation(A, b)






