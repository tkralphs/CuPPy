__version__    = '0.5.2'
__author__     = 'Aykut Bulut and Ted Ralphs'
__license__    = 'Eclipse Public License'
__maintainer__ = 'Ted Ralphs'
__email__      = 'ted@lehigh.edu'
__url__        = 'https://github.com/tkralphs/CuPPy'

import sys
import math
import numpy as np
from cylp.py.modeling import CyLPArray
from milpInstance import MILPInstance

DISPLAY_ENABLED = True
try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    DISPLAY_ENABLED = False

sys.path.append('examples')

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
    f.set_xlim(p.xlim)
    f.set_ylim(p.ylim)
    pI = p.make_integer_hull()
    f.add_polyhedron(pI, show_int_points = True, color = 'red',
                     linestyle = 'dashed',
                     label = 'Convex hull of integer points')
    for (coeff, r) in cuts:
        f.add_line(coeff, r, p.xlim, p.ylim, color = 'green', linestyle = 'dashed')
    if sol is not None:
        f.add_point(sol, radius = .05)
    f.show()


def solve(m, whichCuts = [], 
          debug_print = False, epsilon = .01, 
          max_iter = 100, max_cuts = 10, display = False):    

    if not isinstance(m, MILPInstance):
        raise "Invalid first parameter: Must be of type MILPInstance"

    if not DISPLAY_ENABLED:
        display = False
        
    if m.lp.nCols > 2 or m.A is None:
        display = False
    m.lp.logLevel = 0
    
    if display:
        disp_relaxation(m.A, m.b)
    
    for i in xrange(max_iter):
        print 'Iteration ', i
        m.lp.primal(startFinishOptions = 'x')
        print 'Current bound:', m.lp.objectiveValue
        #Binv = np.zeros(shape = (lp.nConstraints, lp.nConstraints))
        #for i in range(lp.nVariables, lp.nVariables+lp.nConstraints):
        #    lp.getBInvACol(i, Binv[i-lp.nVariables,:])
        #rhs = lp.rhs
        if m.sense == '<=':
            rhs = np.dot(m.lp.basisInverse, m.lp.constraintsUpper)
        else:
            rhs = np.dot(m.lp.basisInverse, m.lp.constraintsLower)
        sol = m.lp.primalVariableSolution['x']
        if debug_print:
            print 'Current basis inverse:'
            print m.lp.basisInverse
            print 'Condition number of basis inverse'
            print np.linalg.cond(m.lp.basisInverse)
            print "Current tableaux:"
            print m.lp.tableau
            print "Current right hand side:\n", rhs
            #print lp.rhs
        print 'Current solution: ', sol
        if isInt(sol[m.integerIndices], epsilon):
            print 'Integer solution found!'
            break
        cuts = []
        for (cg, args) in whichCuts:
            cuts += cg(m.lp, m.integerIndices, m.sense, sol, **args)
        if cuts == []:
            print 'No cuts found!'
            break
        if display:
            disp_relaxation(m.A, m.b, cuts, sol)
        for (coeff, r) in cuts[:max_cuts]:
            #TODO sort cuts by degree of violation
            if m.sense == '<=':
                print 'Adding cut: ', coeff, '<=', r
                m.lp += CyLPArray(coeff) * m.x <= r
            else:
                print 'Adding cut: ', coeff, '>=', r
                m.lp += CyLPArray(coeff) * m.x >= r
            if display:
                m.A.append(coeff.tolist())
                m.b.append(r)
    
    if display:
        disp_relaxation(m.A, m.b)

if __name__ == '__main__':
        
    solve(MILPInstance(module_name = 'MIP6'), whichCuts = [(gomoryCut, {})], display = True, debug_print = True)




