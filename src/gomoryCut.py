import sys
import math
import numpy as np
from CyLP.cy import CyClpSimplex
from CyLP.py.modeling import CyLPArray
from grumpy.polyhedron2D import Polyhedron2D, add_line
import matplotlib.pyplot as plt

epsilon = 0.01
maxiter = 100

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

def gomoryCut(lp, rowInd, value, sense):
    'Return the Gomory cut of row ``rowInd`` of lp (a CyClpSimplex object)'
    f0 = getFraction(value)
    f = np.array([getFraction(x) for x in lp.tableau[rowInd,:lp.nVariables]])
    pi = np.array([f[j]/f0 if f[j] <= f0 
                   else (1-f[j])/(1-f0) for j in range(lp.nVariables)])
    pi_slacks = np.array([x/f0 if x > 0 else -x/(1-f0)  
                         for x in lp.tableau[rowInd,lp.nVariables:]])
    pi -= pi_slacks * lp.coefMatrix
    pi0 = (1 - np.dot(pi_slacks, lp.constraintsUpper) if sense[1] == '<='
           else 1 + np.dot(pi_slacks, lp.constraintsUpper))
    return pi, pi0

def disp_relaxation(A, b):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    p = Polyhedron2D(A = A, b = b)
    p.draw(ax, color = 'blue', linestyle = 'solid')
    ax.set_xlim(p.plot_min[0], p.plot_max[0])
    ax.set_ylim(p.plot_min[1], p.plot_max[1])
    pI = p.make_integer_hull()
    pI.draw(ax, color = 'red', linestyle = 'dashed')
    plt.show()

max_cuts = 1
display = False

lp = CyClpSimplex()

sys.path.append('../instances')

from MIP2 import numVars, numCons, A, b, c, sense, integerIndices
try:
    from MIP2 import x_u
except ImportError:
    x_u = None
else:
    x_u = CyLPArray(x_u)

if display:
    disp_relaxation(A, b)

myA = A
A = np.matrix(A)
#print np.linalg.cond(A)
myb = b
b = CyLPArray(b)

#We assume variables have zero lower bounds
x_l = CyLPArray([0 for i in range(numVars)])

#Integrality tolerance
epsilon = .01

x = lp.addVariable('x', numVars)

lp += x >= x_l
if x_u is not None:
    lp += x <= x_u

lp += (A * x <= b if sense[1] == '<=' else
       A * x >= b)

c = CyLPArray(c)
lp.objective = -c * x if sense[0] == 'Max' else c * x

for i in xrange(maxiter):
    print '______________________________________________________'
    print 'iteration ', i
    lp.primal(startFinishOptions = 'x')
    bv = lp.basicVariables

#    B = np.zeros(shape = (lp.nConstraints, lp.nConstraints))
#    for i in range(lp.nVariables, lp.nVariables+lp.nConstraints):
#        lp.getBInvACol(i, B[i-lp.nVariables,:])
    print 'Current basis inverse:'
    print lp.basisInverse
    print np.linalg.cond(lp.basisInverse)

    print "Current tableaux:"
    print lp.tableau
    print "Current right hand side:"
    print lp.rhs    
    rhs = lp.rhs
    sol = lp.primalVariableSolution['x']
    print 'Current solution: ', sol
    if isInt(sol[integerIndices]):
        print 'Integer solution found!'
        break
    num_cuts = 0
    for rowInd in xrange(lp.nConstraints):
        if num_cuts >= max_cuts:
            break
        basicVarInd = bv[rowInd]
        if basicVarInd in integerIndices and not isInt(sol[basicVarInd]):
            num_cuts += 1
            coef, r = gomoryCut(lp, rowInd, sol[basicVarInd], sense)
            if sense[1] == '>=':
                print 'Adding cut: ', coef, '>=', r
                lp += CyLPArray(coef) * x >= r
                myA.append(coef.tolist())
                myb.append(r)
            else:
                print 'Adding cut: ', -coef, '<=', -r
                lp += CyLPArray(-coef) * x <= -r
                myA.append((-coef).tolist())
                myb.append(-r)
            print myA, myb
            if display:
                disp_relaxation(myA, myb)







