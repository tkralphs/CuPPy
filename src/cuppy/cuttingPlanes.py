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
from math import floor
import numpy as np
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
        return abs(floor(x + 10**(-eps)) - x) < 10**(-eps)
    return (np.abs(np.around(x) - x) < 10**(-eps)).all()

def getFraction(x, eps = EPS):
    'Return the fraction part of x: x - floor(x)'
#    return x - floor(x)
    return np.around(x, decimals = eps) - floor(np.around(x, decimals = eps))

def gomoryMixedIntegerCut(m, rowInds = None, eps = EPS, debug_print = False):
    '''Return the Gomory mixed integer cut of rows in ``rowInds`` of lp 
    (a CyClpSimplex object)'''

    cuts = []
    lp = m.lp
    sol = lp.primalVariableSolution['x']
    if rowInds is None:
        rowInds = list(range(lp.nConstraints))
    for row in rowInds:
        basicVarInd = lp.basicVariables[row]
        if (basicVarInd in m.integerIndices) and (not isInt(sol[basicVarInd], eps)):
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
            pi0 = (1 - np.dot(pi_slacks, lp.constraintsUpper) if m.sense == '<='
                   else 1 + np.dot(pi_slacks, lp.constraintsUpper))
            if m.sense == '>=':
                cuts.append((pi, pi0))
            else:
                cuts.append((-pi, -pi0))
    return cuts, []

def liftAndProject(m, rowInds = None, eps = EPS, debug_print = False):
    '''Return the lift-and-project associated with variables that are basic in 
    rows in ``rowInds`` of lp (a CyClpSimplex object)'''

    cuts = []
    lp = m.lp
    sol = lp.primalVariableSolution['x']
    if rowInds is None:
        rowInds = list(range(lp.nConstraints))
    for row in rowInds:
        basicVarInd = lp.basicVariables[row]
        if (basicVarInd in m.integerIndices) and (not isInt(sol[basicVarInd], eps)): 
            e = np.zeros(lp.nCols)
            e[basicVarInd] = 1
            #Call function for solving CGLP for the associated variable
            #disjunction (disjunction is "<=") 
            cuts += disjunctionToCut(m, e, floor(sol[basicVarInd]), eps = eps,
                                     debug_print = debug_print)
    return cuts, []
            
def disjunctionToCut(m, pi, pi0, debug_print = False, use_cylp = True, eps = EPS):
    '''Generate the most violated valid inequality from a given disjunction'''
    me = "cglp_cuts: "
    lp = m.lp
    sol = lp.primalVariableSolution['x']

    if debug_print:
        print(me, "constraints sense = ", m.sense)
        print(me, "matrix = ")
        print(m.A)
        print(me, "rhs = ", m.b)
        print(me, "vars lower bounds = ", lp.variablesLower)
        print(me, "vars upper bounds = ", lp.variablesUpper)
        print(me, "objective = ", lp.objective)
        print(me, "current solution = ", sol)
        print(me, "pi = ", pi)
        print(me, "pi0 = ", pi0)

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

    pi = CyLPArray(pi)
    
    Atran = m.A.transpose()
    b = CyLPArray(m.b)
    numRows, numCols = m.A.shape
    
    if use_cylp:
        sp = CyLPModel()
        u = sp.addVariable('u', numRows, isInt = False)
        v = sp.addVariable('v', numRows, isInt = False)
        u0 = sp.addVariable('u0', 1, isInt = False)
        v0 = sp.addVariable('v0', 1, isInt = False)
        alpha = sp.addVariable('alpha', lp.nVariables, isInt = False)
        beta = sp.addVariable('beta', 1, isInt = False)
        
        #This should be as simple as this, but it doesn't work.
        #Maybe a bug in CyLP? 
        #sp += alpha - Atran*u - pi*u0 == 0
        #sp += alpha - Atran*v + pi*v0 == 0
        for i in range(numCols):
            sp += alpha[i] - sum(Atran[i,j]*u[j] for j in range(numRows)) - pi[i]*u0 == 0
        for i in range(numCols):
            sp += alpha[i] - sum(Atran[i,j]*v[j] for j in range(numRows)) + pi[i]*v0 == 0
        if m.sense == '<=':
            sp += beta - b*u - pi0*u0 >= 0
            sp += beta - b*v + (pi0 + 1)*v0 >= 0
        else:
            sp += beta - b*u - pi0*u0 <= 0
            sp += beta - b*v + (pi0 + 1)*v0 <= 0
        sp += u0 + v0 == 1
        sp += u >= 0
        sp += v >= 0
        sp += u0 >= 0
        sp += v0 >= 0
        if m.sense == '<=':
            sp.objective = sum(-sol[i]*alpha[i] for i in range(numCols)) + beta
        else:
            #This direction is not debugged
            sp.objective = sum(sol[i]*alpha[i] for i in range(numCols)) - beta            

        cglp = CyClpSimplex(sp)
        # If we want to solve it as an MILP
        # cglp = CyClpSimplex(sp).getCbcModel()
        #cglp.writeLp('lp.lp')
        cglp.logLevel = 0
        cglp.primal(startFinishOptions = 'x')
        # Solve as MILP
        # cglp.solve()
        beta = cglp.primalVariableSolution['beta'][0]
        alpha = cglp.primalVariableSolution['alpha']
        u = cglp.primalVariableSolution['u']
        v = cglp.primalVariableSolution['v']
        u0 = cglp.primalVariableSolution['u0'][0]
        v0 = cglp.primalVariableSolution['v0'][0]
        if debug_print:
            print(me, 'Objective Value: ', cglp.objectiveValue)

        if debug_print:
            print(me, 'u: ', u)
            print(me, 'v: ', v)
            print(me, 'u0: ', u0)
            print(me, 'v0: ', v0)
    else: 
        CG = AbstractModel()
        CG.u = Var(list(range(numRows)), domain=NonNegativeReals,
                   bounds = (0.0, None))
        CG.v = Var(list(range(numRows)), domain=NonNegativeReals,
                   bounds = (0.0, None))
        CG.u0 = Var(domain=NonNegativeReals, bounds = (0.0, None))
        CG.v0 = Var(domain=NonNegativeReals, bounds = (0.0, None))
        CG.alpha = Var(list(range(numRows)), domain=Reals,
                       bounds = (None, None))    
        CG.beta  = Var(domain=Reals, bounds = (None, None))    
        
        ## Constraints
        def pi_rule_left(CG, i):
            x = float(pi[i])
            return(sum(Atran[i, j]*CG.u[j] for j in range(numRows)) -
                   x*CG.u0 - CG.alpha[i] == 0.0)
        CG.pi_rule_left = Constraint(list(range(numCols)), rule=pi_rule_left)
        
        def pi_rule_right(CG, i):
            x = float(pi[i])
            return(sum(Atran[i, j]*CG.v[j] for j in range(numRows)) +
                   x*CG.v0 - CG.alpha[i] == 0.0)
        CG.pi_rule_right = Constraint(list(range(numCols)), rule=pi_rule_right)

        if m.sense == '<=':
            def pi0_rule_left(CG):
                return(sum(b[j]*CG.u[j] for j in range(numRows)) -
                       pi0*CG.u0 - CG.beta <= 0.0)
            CG.pi0_rule_left = Constraint(rule=pi0_rule_left)
            
            def pi0_rule_right(CG):
                return(sum(b[j]*CG.v[j] for j in range(numRows)) +
                       (pi0 + 1)*CG.v0 - CG.beta <= 0.0)
            CG.pi0_rule_right = Constraint(rule=pi0_rule_right)
        else:
            def pi0_rule_left(CG):
                return(sum(b[j]*CG.u[j] for j in range(numRows)) -
                       pi0*CG.u0 - CG.beta >= 0.0)
            CG.pi0_rule_left = Constraint(rule=pi0_rule_left)
            
            def pi0_rule_right(CG):
                return(sum(b[j]*CG.v[j] for j in range(numRows)) +
                       (pi0 + 1)*CG.v0 - CG.beta >= 0.0)
            CG.pi0_rule_right = Constraint(rule=pi0_rule_right)

        def normalization_rule(CG):
            return(CG.u0 + CG.v0 == 1.0)
        CG.normalization_rule = Constraint(rule=normalization_rule)
        
        def objective_rule(CG):
            return(sum(sol[i]*CG.alpha[i] for i in range(numCols)) -
                   CG.beta)
        if m.sense == '<=':
            CG.objective = Objective(sense=maximize, rule=objective_rule)
        else:
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

    
    if np.abs(violation) > 10**(-eps):
        return [(alpha, beta)]

    print('No violated cuts found solving CGLP', violation)
    return []

def disp_relaxation(f, A, b, cuts = [], sol = None, disj = [], filename = None):
    #TODO: Check sense of inequalities by looking explicitly at
    #      lp.constraintsUpper and lp.constraintsLower
    p = Polyhedron2D(A = A, b = b)
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
    f.show(filename = filename)

def solve(m, whichCuts = [], use_cglp = False, debug_print = False, eps = EPS, 
          max_iter = 100, max_cuts = 10, display = False, filename = None):    

    if not isinstance(m, MILPInstance):
        print("Invalid first parameter: Must be of type MILPInstance")
        exit

    if not DISPLAY_ENABLED:
        display = False
    else:
        f = Figure()
  
    if m.lp.nCols > 2 or m.A is None:
        display = False
    m.lp.logLevel = 0

    #Include bounds explicitly in the constraint matrix for display and for
    #use in cut generators. 
    infinity = m.lp.getCoinInfinity()
    if m.sense == '<=':
        b = m.lp.constraintsUpper.copy()
        mult = -1.0
    else:
        b = m.lp.constraintsLower.copy()
        mult = 1.0
    if type(m.A) == csc_matrixPlus:
        A = m.A.toarray()
    else:
        A = m.A.copy()
    for i in range(m.lp.nCols):
        e = np.zeros((1, m.lp.nCols))
        if m.lp.variablesUpper[i] < infinity:
            b.resize(b.size+1, refcheck = False)
            e[0, i] = -mult
            b[-1] = -mult*m.lp.variablesUpper[i]
            A = np.vstack((A, e))
        if m.lp.variablesLower[i] > -infinity:
            b.resize(b.size+1, refcheck = False)
            e[0, i] = mult
            b[-1] = mult*m.lp.variablesLower[i]
            A = np.vstack((A, e))
    m.A = A
    m.b = b

    if display and filename is not None:
        disp_relaxation(f, m.A, m.b, filename = filename+'.png')
    elif display:
        disp_relaxation(f, m.A, m.b)

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
            print('Condition number of basis inverse',
                  np.around(np.linalg.cond(m.lp.basisInverse)))
            print('Current tableaux:')
            print(m.lp.tableau)
            print('Current right hand side:\n', rhs)
            #print('Dual solution:', m.lp.dualConstraintSolution)
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
                tmp_cuts, tmp_disj = cg(m, **args, eps = eps)
                cuts += tmp_cuts
                disj += tmp_disj
        cur_num_cuts = len(cuts)
        if use_cglp:
            if len(disj) > 0:
                for d in disj:
                    cuts += disjunctionToCut(m, d[0], d[1], eps = eps)
        if cuts == []:
            if disj == []:
                print('No cuts found and terminating!')
                break
            else:
                print('No cuts found but continuing!')
        if display and filename is not None:
            disp_relaxation(f, m.A, m.b, cuts, sol, disj,
                            filename = filename+str(i)+'.png')
        elif display:
            disp_relaxation(f, m.A, m.b, cuts, sol, disj)
        if len(cuts) == cur_num_cuts:
            disj = []
        for (coeff, r) in cuts[:max_cuts]:
            #TODO sort cuts by degree of violation
            if m.sense == '<=':
                coeff = np.floor(coeff*(10**eps))/(10**eps)
                r = np.ceil(r*(10**eps))/(10**eps)
                print('Adding cut: ', coeff, '<=', r)
                m.lp += CyLPArray(coeff) * m.x <= r
            else:
                coeff = np.ceil(coeff*(10**eps))/(10**eps)
                r = np.floor(r*(10**eps))/(10**eps)
                print('Adding cut: ', coeff, '>=', r)
                m.lp += CyLPArray(coeff) * m.x >= r
            m.A = np.vstack((m.A, np.array(coeff)))
            m.b.resize(m.b.size+1, refcheck = False)
            m.b[-1] = r
            
    if display:
        disp_relaxation(f, m.A, m.b)

if __name__ == '__main__':
            
    solve(MILPInstance(module_name = 'coinor.cuppy.examples.MIP6'),
          whichCuts = [(gomoryMixedIntegerCut, {})],
          display = True, debug_print = True, use_cglp = False)




