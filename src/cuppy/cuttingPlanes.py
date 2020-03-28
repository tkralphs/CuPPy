from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from past.utils import old_div
__version__    = '0.5.2'
__author__     = 'Aykut Bulut and Ted Ralphs'
__license__    = 'Eclipse Public License'
__maintainer__ = 'Ted Ralphs'
__email__      = 'ted@lehigh.edu'
__url__        = 'https://github.com/tkralphs/CuPPy'

import sys
import math
import numpy as np
from copy import deepcopy
from cylp.py.utils.sparseUtil import csc_matrixPlus 
from cylp.cy import CyClpSimplex
from cylp.py.modeling import CyLPArray, CyLPModel
PYOMO_INSTALLED = True
try:
    from pyomo.environ import AbstractModel, Var, Constraint, SolverFactory
    from pyomo.environ import NonNegativeReals, NonPositiveReals, Reals, Set
    from pyomo.environ import Integers, Objective, minimize, value
except ImportError:
    PYOMO_INSTALLED = False

DISPLAY_ENABLED = True
try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
    from src.cuppy.milpInstance import MILPInstance
except ImportError:
    try:
        from coinor.grumpy.polyhedron2D import Polyhedron2D, Figure
        from coinor.cuppy.milpInstance import MILPInstance
    except ImportError:
        DISPLAY_ENABLED = False

sys.path.append('examples')

EPS = 5

def isInt(x, eps = EPS):
    '''
    Return True if x is an integer, or if x is a numpy array
    with all integer elements, False otherwise
    '''
    if isinstance(x, (int, float)):
        return abs(math.floor(x + 10**(-eps)) - x) < 10**(-eps)
    return (np.abs(np.around(x) - x) < 10**(-eps)).all()

def getFraction(x, eps = EPS):
    'Return the fraction part of x: x - floor(x)'
#    return x - math.floor(x)
    return np.around(x, decimals = eps) - math.floor(np.around(x, decimals = eps))

def gomoryCut(lp, integerIndices = None, sense = '>=', sol = None,
              rowInds = None, value = None, eps = EPS):
    '''Return the Gomory cut of rows in ``rowInds`` of lp 
    (a CyClpSimplex object)'''
    cuts = []
    if sol is None:
        sol = lp.primalVariableSolution['x']
    if rowInds is None:
        rowInds = list(range(lp.nConstraints))
    if integerIndices is None:
        integerIndices = list(range(lp.nVariables))
    for row in rowInds:
        basicVarInd = lp.basicVariables[row]
        if (basicVarInd in integerIndices) and (not isInt(sol[basicVarInd], eps)):
            f0 = getFraction(sol[basicVarInd], eps)
            f = []
            for i in range(lp.nVariables):
                if i in lp.basicVariables:
                    #This is to try to avoid getting very small numbers that 
                    #should be zero
                    f.append(0)
                else:
                    f.append(getFraction(lp.tableau[row, i], eps))
            pi = np.array([old_div(f[j],f0) if f[j] <= f0 
                           else old_div((1-f[j]),(1-f0)) for j in range(lp.nVariables)])
            pi_slacks = np.array([old_div(x,f0) if x > 0 else old_div(-x,(1-f0))  
                                 for x in lp.tableau[row, lp.nVariables:]])
            pi -= pi_slacks * lp.coefMatrix
            pi0 = (1 - np.dot(pi_slacks, lp.constraintsUpper) if sense == '<='
                   else 1 + np.dot(pi_slacks, lp.constraintsUpper))
            pi = np.ceil(pi*(10**eps))/(10**eps)
            pi0 = np.floor(pi0*(10**eps))/(10**eps)
            if sense == '>=':
                cuts.append((pi, pi0))
            else:
                cuts.append((-pi, -pi0))
    return cuts, []
            
def disjunctionToCut(lp, pi, pi0, integerIndices = None, sense = '>=',
                     sol = None, debug_print = False, use_cylp = True,
                     eps = EPS):

    me = "cglp_cuts: "

    if sol is None:
        sol = lp.primalVariableSolution['x']
    infinity = lp.getCoinInfinity()

    if debug_print:
        print(me, "constraints sense = ", sense)
        print(me, "con lower bounds = ", lp.constraintsLower)
        print(me, "con upper bounds = ", lp.constraintsUpper)
        print(me, "con matrix = ", lp.coefMatrix.toarray())
        print(me, "vars lower bounds = ", lp.variablesLower)
        print(me, "vars upper bounds = ", lp.variablesUpper)
        print(me, "Assuming objective is to minimize")
        print(me, "objective = ", lp.objective)
        print(me, "infinity = ", infinity)
        print(me, "current point = ", sol)
        print(me, "pi = ", pi)
        print(me, "pi0 = ", pi0)

    A = lp.coefMatrix.toarray()
    #c = lp.objective
    ## Convert to >= if the problem is in <= form.
    if sense == '<=':
        b = deepcopy(lp.constraintsUpper)
        b = -1.0*b
        A = -1.0*A
    else:
        b = deepcopy(lp.constraintsLower)

    #Add bounds on variables as explicit constraints
    for i in range(lp.nCols):
        e = np.zeros((1, lp.nCols))
        if lp.variablesUpper[i] < infinity:
            b.resize(b.size+1, refcheck = False)
            e[0, i] = -1.0
            b[-1] = -1.0*lp.variablesUpper[i]
            A = np.vstack((A, e))
        if lp.variablesLower[i] > -infinity:
            b.resize(b.size+1, refcheck = False)
            e[0, i] = 1.0
            b[-1] = lp.variablesLower[i]
            A = np.vstack((A, e))
    A = csc_matrixPlus(A)

    ############################################################################
    ## There are two given LPs:
    ## s.t. Ax >= b           s.t. Ax >= b
    ##   -pi.x >= -pi_0          pi.x >= pi_0+1
    ## A, b, c, pi, pi_0 are given
    ##
    ## CGLP: alpha.x >= beta should be valid for both LPs above
    ##
    ## min alpha.x* - beta
    ## uA - u0.pi = alpha
    ## vA + v0.pi = alpha
    ## ub - u0.pi_0 >= beta 
    ## vb + v0.(pi_0 + 1) >= beta 
    ## u0 + v0 = 1
    ## u, v, u0, v0 >= 0
    ## if min value comes out < 0, then (alpha.x >= beta) is a cut.
    ############################################################################

    b = CyLPArray(b)
    pi = CyLPArray(pi)
    
    Atran = A.transpose()

    if use_cylp:
        sp = CyLPModel()
        u = sp.addVariable('u', A.shape[0], isInt = False)
        v = sp.addVariable('v', A.shape[0], isInt = False)
        u0 = sp.addVariable('u0', 1, isInt = False)
        v0 = sp.addVariable('v0', 1, isInt = False)
        alpha = sp.addVariable('alpha', lp.nVariables, isInt = False)
        beta = sp.addVariable('beta', 1, isInt = False)
    
        for i in range(A.shape[1]):
            sp += alpha[i] - sum(Atran[i,j]*u[j] for j in range(A.shape[0])) + pi[i]*u0 == 0
        for i in range(A.shape[1]):
            sp += alpha[i] - sum(Atran[i,j]*v[j] for j in range(A.shape[0])) - pi[i]*v0 == 0
        sp += beta - b*u + pi0*u0 <= 0
        sp += beta - b*v - (pi0 + 1)*v0 <= 0
        sp += u0 + v0 == 1
        if sense == '<=':
            sp += u >= 0
            sp += v >= 0
            sp += u0 >= 0
            sp += v0 >= 0
        else:
            #TODO this direction is not debugged
            # Is this all we need?
            sp += u <= 0
            sp += v <= 0
            sp += u0 <= 0
            sp += v0 <= 0
        sp.objective = sum(sol[i]*alpha[i] for i in range(A.shape[1])) - beta
        cbcModel = CyClpSimplex(sp).getCbcModel()
        cbcModel.logLevel = 0
        #cbcModel.maximumSeconds = 5
        cbcModel.solve()
        beta = cbcModel.primalVariableSolution['beta'][0]
        alpha = cbcModel.primalVariableSolution['alpha']
        u = cbcModel.primalVariableSolution['u']
        v = cbcModel.primalVariableSolution['v']
        u0 = cbcModel.primalVariableSolution['u0'][0]
        v0 = cbcModel.primalVariableSolution['v0'][0]
        if debug_print:
            print(me, 'Objective Value: ', cbcModel.objectiveValue)
            print(me, 'alpha: ', alpha, 'alpha*sol: ', np.dot(alpha, sol))
            print(me, 'beta: ', beta)
            print(me, 'Violation of cut: ',  np.dot(alpha, sol) - beta)
    else: 
        CG = AbstractModel()
        CG.u = Var(list(range(A.shape[0])), domain=NonNegativeReals,
                   bounds = (0.0, None))
        CG.v = Var(list(range(A.shape[0])), domain=NonNegativeReals,
                   bounds = (0.0, None))
        CG.u0 = Var(domain=NonNegativeReals, bounds = (0.0, None))
        CG.v0 = Var(domain=NonNegativeReals, bounds = (0.0, None))
        CG.alpha = Var(list(range(A.shape[0])), domain=Reals,
                       bounds = (None, None))    
        CG.beta  = Var(domain=Reals, bounds = (None, None))    
        
        ## Constraints
        def pi_rule_left(CG, i):
            x = float(pi[i])
            return(sum(Atran[i, j]*CG.u[j] for j in range(A.shape[0])) -
                   x*CG.u0 - CG.alpha[i] == 0.0)
        CG.pi_rule_left = Constraint(list(range(A.shape[1])), rule=pi_rule_left)
        
        def pi_rule_right(CG, i):
            x = float(pi[i])
            return(sum(Atran[i, j]*CG.v[j] for j in range(A.shape[0])) +
                   x*CG.v0 - CG.alpha[i] == 0.0)
        CG.pi_rule_right = Constraint(list(range(A.shape[1])), rule=pi_rule_right)
        
        def pi0_rule_left(CG):
            return(sum(b[j]*CG.u[j] for j in range(A.shape[0])) -
                   pi0*CG.u0 - CG.beta >= 0.0)
        CG.pi0_rule_left = Constraint(rule=pi0_rule_left)
        
        def pi0_rule_right(CG):
            return(sum(b[j]*CG.v[j] for j in range(A.shape[0])) +
                   (pi0 + 1)*CG.v0 - CG.beta >= 0.0)
        CG.pi0_rule_right = Constraint(rule=pi0_rule_right)
        
        def normalization_rule(CG):
            return(CG.u0 + CG.v0 == 1.0)
        CG.normalization_rule = Constraint(rule=normalization_rule)
        
        def objective_rule(CG):
            return(sum(sol[i]*CG.alpha[i] for i in range(A.shape[1])) -
                   CG.beta)
        CG.objective = Objective(sense=minimize, rule=objective_rule)
        
        opt = SolverFactory("cbc")
        instance = CG.create_instance()
        #instance.pprint()
        #instance.write("foo.nl", format = "nl")
        #opt.options['bonmin.bb_log_level'] = 5
        #opt.options['bonmin.bb_log_interval'] = 1
        results = opt.solve(instance, tee=False)
        #results = opt.solve(instance)
        instance.solutions.load_from(results)
        
        beta = instance.beta.value
        alpha = np.array([instance.alpha[i].value
                          for i in range(lp.nVariables)])
    violation =  beta - np.dot(alpha, sol) 
    if debug_print:
        print(me, 'Beta: ', beta)
        print(me, 'alpha: ', alpha)
        print(me, 'Violation of cut: ', violation)
        
    if violation > 10**(-eps):
        if (sense == ">="):
            return [(alpha, beta)]
        else:
            return [(-alpha, -beta)]
    return []

def disp_relaxation(A, b, cuts = [], sol = None, disj = []):
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
    for (coeff, r) in disj:
        f.add_line(coeff, r, p.xlim, p.ylim, color = 'red', linestyle = 'dashed')
        f.add_line(coeff, r+1, p.xlim, p.ylim, color = 'red', linestyle = 'dashed')
    if sol is not None:
        f.add_point(sol, radius = .05)
    f.show()


def solve(m, whichCuts = [], use_cglp = False,
          debug_print = False, eps = EPS, 
          max_iter = 100, max_cuts = 10, display = False):    

    if not isinstance(m, MILPInstance):
        print("Invalid first parameter: Must be of type MILPInstance")
        exit

    if not DISPLAY_ENABLED:
        display = False
        
    if m.lp.nCols > 2 or m.A is None:
        display = False
    m.lp.logLevel = 0
    
    if display:
        disp_relaxation(m.A, m.b)
    
    disj = []
    prev_sol = np.zeros((1, m.lp.nCols))
    for i in range(max_iter):
        print('Iteration ', i)
        m.lp.primal(startFinishOptions = 'x')
        print('Current bound:', m.lp.objectiveValue)
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
            print('Current basis inverse:')
            print(m.lp.basisInverse)
            print('Condition number of basis inverse', np.around(np.linalg.cond(m.lp.basisInverse)))
            print('Current tableaux:')
            print(m.lp.tableau)
            print('Current right hand side:\n', rhs)
            #print lp.rhs
        print('Current solution: ', sol)

        if (sol - prev_sol).any():
            prev_sol = sol
        else:
            print ("Solution repeated, stalling detected") 
            print ("Exiting")
            break

        if isInt(sol[m.integerIndices], eps):
            print('Integer solution found!')
            break

        if np.around(np.linalg.cond(m.lp.basisInverse)) >= 10**32:
            print ("Condition number of the basis matrix exceeds 10^32") 
            print ("Exiting")
            break

        cuts = []
        if disj == []:
            for (cg, args) in whichCuts:
                tmp_cuts, tmp_disj = cg(m.lp, m.integerIndices, m.sense, sol, **args, eps = eps)
                cuts += tmp_cuts
                disj += tmp_disj
        cur_num_cuts = len(cuts)
        if use_cglp and len(disj) > 0:
            for d in disj:
                cuts += disjunctionToCut(m.lp, d[0], d[1], sense=m.sense, eps = eps)
        if cuts == []:
            if disj == []:
                print('No cuts found and terminating!')
                break
            else:
                print('No cuts found but continuing!')
        if display:
            disp_relaxation(m.A, m.b, cuts, sol, disj)
        if len(cuts) == cur_num_cuts:
            disj = []
        for (coeff, r) in cuts[:max_cuts]:
            #TODO sort cuts by degree of violation
            if m.sense == '<=':
                print('Adding cut: ', coeff, '<=', r)
                m.lp += CyLPArray(coeff) * m.x <= r
            else:
                print('Adding cut: ', coeff, '>=', r)
                m.lp += CyLPArray(coeff) * m.x >= r
            if display:
                m.A.append(coeff.tolist())
                m.b.append(r)
    
    if display:
        disp_relaxation(m.A, m.b)

if __name__ == '__main__':
            
    solve(MILPInstance(module_name = 'coinor.cuppy.examples.MIP6'), whichCuts = [(gomoryCut, {})], display = True, debug_print = True,
          use_cglp = False)




