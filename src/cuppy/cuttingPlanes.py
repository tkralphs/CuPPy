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
PYOMO_INSTALLED = True
try:
    from pyomo.environ import *
except ImportError:
    PYOMO_INSTALLED = False

sys.path.append('instances')

DISPLAY_ENABLED = True
try:
    from src.grumpy.polyhedron2D import Polyhedron2D, Figure
except ImportError:
    DISPLAY_ENABLED = False

def isInt(x):
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

def maxViolationSplitCuts(lp, integerIndices = None, sense = '>=', sol = None,
              max_coeff = 1):
    #Warning: At the moment, you must put bound constraints in explicitly for split cuts
    A = lp.coefMatrix
    if sense == '<=':
        b = CyLPArray(lp.constraintsUpper)
    else:
        b = CyLPArray(lp.constraintsLower)
    if integerIndices is None:
        integerIndices = range(lp.nVariables)
    if sol is None:
        sol = lp.primalVariableSolution['x']
    s = A*sol - b
    best = lp.getCoinInfinity()
    best_theta = None

    for theta in [0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5]:
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
        cbcModel.logLevel = 0
        #cbcModel.maximumSeconds = 5
        cbcModel.solve()
        if debug_print:
            #print 'Theta: ', theta, 
            #print 'Objective Value: ', cbcModel.objectiveValue - theta*(1-theta)
            #print 'pi: ', cbcModel.primalVariableSolution['pi']
            #print 'pi0: ', cbcModel.primalVariableSolution['pi0']
            multu = cbcModel.primalVariableSolution['u']
            disjunction = cbcModel.primalVariableSolution['pi']
            rhs = cbcModel.primalVariableSolution['pi0']
            alpha = A.transpose()*multu + theta*disjunction
            beta = np.dot(b, multu) + theta*rhs
            #print 'alpha: ', alpha, 'alpha*sol: ', np.dot(alpha, sol)
            #print 'beta: ', beta
            #print 'Violation of cut: ',  np.dot(alpha, sol) - beta
        if cbcModel.objectiveValue - theta*(1-theta) < best:
            best = cbcModel.objectiveValue - theta*(1-theta)
            best_multu = cbcModel.primalVariableSolution['u']
            best_multv = cbcModel.primalVariableSolution['v']
            best_disjunction = cbcModel.primalVariableSolution['pi']
            best_rhs = cbcModel.primalVariableSolution['pi0']
            best_theta = theta
            
    if best_theta is not None:
        alpha = A.transpose()*best_multu + best_theta*best_disjunction
        beta = np.dot(b, best_multu) + best_theta*best_rhs
        if debug_print:
            print 'Violation of cut: ',  np.dot(alpha, sol) - beta
            print 'pi: ', best_disjunction
            print 'pi0: ', best_rhs
            print 'theta: ',  best_theta
        if (abs(alpha) > 1e-6).any():
            return [(alpha, beta)] 
    return []

def maxViolationSplitCuts2(lp, integerIndices = None, sense = '>=', sol = None,
              max_coeff = 1):
    #Warning: At the moment, you must put bound constraints in explicitly for split cuts
    A = lp.coefMatrix
    if sense == '<=':
        b = CyLPArray(lp.constraintsUpper)
    else:
        b = CyLPArray(lp.constraintsLower)
    if integerIndices is None:
        integerIndices = range(lp.nVariables)
    if sol is None:
        sol = lp.primalVariableSolution['x']
    s = A*sol - b

    CG = AbstractModel()
    CG.intIndices = Set(initialize=integerIndices)
    if sense == '<=':
        CG.u = Var(range(lp.nConstraints), domain=NonNegativeReals, bounds = (None, None))
        CG.v = Var(range(lp.nConstraints), domain=NonNegativeReals, bounds = (None, None))
    else:        
        #TODO this direction is not debugged
        # Is this all we need?
        CG.u = Var(range(lp.nConstraints), domain=NonPositiveReals, bounds = (None, None))
        CG.v = Var(range(lp.nConstraints), domain=NonPositiveReals, bounds = (None, None))
    #This assumes that all variables are integer valued for now
    CG.pi = Var(range(lp.nVariables), domain=Integers, bounds = (-max_coeff, max_coeff))
    CG.pi0 = Var([1], domain=Integers, bounds = (None, None))    
    CG.theta = Var([1], domain=NonNegativeReals, bounds = (0, 0.5))
    Atran = A.transpose()

    def pi_rule(CG, i):
        return(CG.pi[i] + sum(Atran[i, j]*CG.u[j] - Atran[i, j]*CG.v[j] for j in range(lp.nConstraints)) == 0)
    CG.pi_rule = Constraint(range(lp.nVariables), rule=pi_rule)
    
    def pi0_rule(CG):
        return(CG.pi0[1] + sum(b[i]*CG.u[i] - b[i]*CG.v[i] for i in range(lp.nConstraints)) ==  CG.theta[1] -1)
        #return(b[0]*CG.v[0] + CG.theta == CG.pi0)
    CG.pi0_rule = Constraint(rule=pi0_rule)
            
    def objective_rule(CG):
        return((CG.theta[1]-1)*sum(s[i]*CG.u[i] for i in range(lp.nConstraints)) - 
               CG.theta[1]*sum(s[i]*CG.v[i] for i in range(lp.nConstraints)) - 
               CG.theta[1]*(1-CG.theta[1]))
    CG.objective = Objective(sense=minimize, rule=objective_rule)

    debug_print = True
    opt = SolverFactory("couenne")
    #opt.options['display_stats'] = 'no' 
    instance = CG.create()
    #instance.pprint()
    #opt.options['bonmin.bb_log_level'] = 5
    #opt.options['bonmin.bb_log_interval'] = 1
    results = opt.solve(instance)
    instance.load(results)

    multu = np.array([instance.u[i].value for i in range(lp.nConstraints)])
    disjunction = np.array([instance.pi[i].value for i in range(lp.nVariables)])
    rhs = instance.pi0[1].value
    theta = instance.theta[1].value
    alpha = A.transpose()*multu + theta*disjunction
    beta = np.dot(b, multu) + theta*rhs
    if debug_print:
        print 'Theta: ', theta, 
        print 'Objective Value: ', value(instance.objective)
        print 'pi: ', disjunction
        print 'pi0: ', rhs
        print 'Violation of cut: ',  np.dot(alpha, sol) - beta

    if (abs(alpha) > 1e-6).any():
        return [(alpha, beta)] 
    return []

def maxViolationSplitCuts3(lp, integerIndices = None, sense = '>=', sol = None,
                 max_coeff = 1):
    #Warning: At the moment, you must put bound constraints in explicitly for split cuts
    A = lp.coefMatrix
    if sense == '<=':
        b = CyLPArray(lp.constraintsUpper)
    else:
        b = CyLPArray(lp.constraintsLower)
    if integerIndices is None:
        integerIndices = range(lp.nVariables)
    if sol is None:
        sol = lp.primalVariableSolution['x']
    s = A*sol - b

    CG = AbstractModel()
    CG.intIndices = Set(initialize=integerIndices)
    if sense == '<=':
        CG.u = Var(range(lp.nConstraints), domain=NonNegativeReals, bounds = (None, None))
        CG.v = Var(range(lp.nConstraints), domain=NonNegativeReals, bounds = (None, None))
    else:        
        #TODO this direction is not debugged
        # Is this all we need?
        CG.u = Var(range(lp.nConstraints), domain=NonPositiveReals, bounds = (None, None))
        CG.v = Var(range(lp.nConstraints), domain=NonPositiveReals, bounds = (None, None))
    #This assumes that all variables are integer valued for now
    CG.pi = Var(range(lp.nVariables), domain=Integers, bounds = (-max_coeff, max_coeff))
    CG.pi0 = Var([1], domain=Integers, bounds = (None, None))    
    CG.theta = 0.5
    Atran = A.transpose()

    def pi_rule(CG, i):
        return(CG.pi[i] + sum(Atran[i, j]*CG.u[j] - Atran[i, j]*CG.v[j] for j in range(lp.nConstraints)) == 0)
    CG.pi_rule = Constraint(range(lp.nVariables), rule=pi_rule)
    
    def pi0_rule(CG):
        return(CG.pi0[1] + sum(b[i]*CG.u[i] - b[i]*CG.v[i] for i in range(lp.nConstraints)) ==  CG.theta -1)
        #return(b[0]*CG.v[0] + CG.theta == CG.pi0)
    CG.pi0_rule = Constraint(rule=pi0_rule)
            
    def objective_rule(CG):
        return((CG.theta-1)*sum(s[i]*CG.u[i] for i in range(lp.nConstraints)) - 
               CG.theta*sum(s[i]*CG.v[i] for i in range(lp.nConstraints)) - 
               CG.theta*(1-CG.theta))
    CG.objective = Objective(sense=minimize, rule=objective_rule)

    debug_print = True
    opt = SolverFactory("couenne")
    #opt.options['bonmin.bb_log_level'] = 5
    #opt.options['bonmin.bb_log_interval'] = 1
    instance = CG.create()
    #instance.write("foo.nl", format="nl")
    #instance.pprint()
    #results = opt.solve(instance, tee=True)
    results = opt.solve(instance)
    instance.load(results)

    multu = np.array([instance.u[i].value for i in range(lp.nConstraints)])
    disjunction = np.array([instance.pi[i].value for i in range(lp.nVariables)])
    rhs = instance.pi0[1].value
    alpha = A.transpose()*multu + CG.theta*disjunction
    beta = np.dot(b, multu) + CG.theta*rhs
    if debug_print:
        print 'Theta: ', CG.theta, 
        print 'Objective Value: ', value(instance.objective)
        print 'pi: ', disjunction
        print 'pi0: ', rhs
        print 'Violation of cut: ',  np.dot(alpha, sol) - beta

    if (abs(alpha) > 1e-6).any():
        return [(alpha, beta)] 
    return []

def boundOptimalSplitCuts(lp, integerIndices = None, sense = '>=', sol = None, 
               max_coeff = 1):
    #Warning: At the moment, you must put bound constraints in explicitly for split cuts
    #This will probably only work with non-negative variables
    A = lp.coefMatrix
    if sense == '<=':
        b = CyLPArray(lp.constraintsUpper)
    else:
        b = CyLPArray(lp.constraintsLower)
    c = lp.objective
    if integerIndices is None:
        integerIndices = range(lp.nVariables)
    if sol is None:
        sol = lp.primalVariableSolution['x']

    CG = AbstractModel()
    CG.intIndices = Set(initialize=integerIndices)
    if sense == '<=':
        CG.u = Var(range(lp.nConstraints), domain=NonPositiveReals, bounds = (None, None))
        CG.v = Var(range(lp.nConstraints), domain=NonPositiveReals, bounds = (None, None))
        CG.u0 = Var([1], domain=NonPositiveReals, bounds = (None, None))
        CG.v0 = Var([1], domain=NonPositiveReals, bounds = (None, None))
    else:        
        #TODO this direction is not debugged
        # Is this all we need?
        CG.u = Var(range(lp.nConstraints), domain=NonNegativeReals, bounds = (None, None))
        CG.v = Var(range(lp.nConstraints), domain=NonNegativeReals, bounds = (None, None))
        CG.u0 = Var([1], domain=NonNegativeReals, bounds = (None, None))
        CG.v0 = Var([1], domain=NonNegativeReals, bounds = (None, None))
    #This assumes that all variables are integer valued for now
    CG.pi = Var(range(lp.nVariables), domain=Integers, bounds = (-max_coeff, max_coeff))
    CG.pi0 = Var([1], domain=Integers, bounds = (None, None))    
    CG.beta = Var([1], domain=Reals, bounds = (None, None))    
    Atran = A.transpose()

    def pi_rule_left(CG, i):
        return(sum(Atran[i, j]*CG.u[j] for j in range(lp.nConstraints))  + CG.u0[1]*CG.pi[i] <= c[i])
    CG.pi_rule_left = Constraint(range(lp.nVariables), rule=pi_rule_left)

    def pi_rule_right(CG, i):
        return(sum(Atran[i, j]*CG.v[j] for j in range(lp.nConstraints)) - CG.v0[1]*CG.pi[i] <= c[i])
    CG.pi_rule_right = Constraint(range(lp.nVariables), rule=pi_rule_right)
    
    def pi0_rule_left(CG):
        return(sum(b[j]*CG.u[j] for j in range(lp.nConstraints)) + CG.u0[1]*CG.pi0[1] >= CG.beta[1])
        #return(b[0]*CG.v[0] + CG.theta == CG.pi0)
    CG.pi0_rule_left = Constraint(rule=pi0_rule_left)

    def pi0_rule_right(CG):
        return(sum(b[j]*CG.v[j] for j in range(lp.nConstraints)) - CG.v0[1]*(CG.pi0[1] + 1) >= CG.beta[1])
        #return(b[0]*CG.v[0] + CG.theta == CG.pi0)
    CG.pi0_rule_right = Constraint(rule=pi0_rule_right)
        
    def objective_rule(CG):
        return(-CG.beta[1])
    CG.objective = Objective(sense=minimize, rule=objective_rule)

    debug_print = True
    opt = SolverFactory("couenne")
    instance = CG.create()
    #instance.pprint()
    #instance.write("foo.nl", format = "nl")
    #opt.options['bonmin.bb_log_level'] = 5
    #opt.options['bonmin.bb_log_interval'] = 1
    #results = opt.solve(instance, tee=True)
    results = opt.solve(instance)
    instance.load(results)

    #multu = np.array([instance.u[i].value for i in range(lp.nConstraints)])
    disjunction = np.array([instance.pi[i].value for i in range(lp.nVariables)])
    rhs = instance.pi0[1].value
    beta = instance.beta[1].value
    violation =  beta - np.dot(c, sol)
    if debug_print:
        print 'Beta: ', beta
        print 'pi: ', disjunction
        print 'pi0: ', rhs
        print 'Violation of cut: ', violation

    if (violation > 1e-6):
        return [(c, beta)] 
    return []
    
def gomoryCut(lp, integerIndices = None, sense = '>=', rowInds = None, 
              value = None):
    '''Return the Gomory cut of rows in ``rowInds`` of lp 
    (a CyClpSimplex object)'''
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
            pass
        
        lp += (A * x <= b if mip.sense[1] == '<=' else
               A * x >= b)
        c = CyLPArray(mip.c)
        lp.objective = -c * x if mip.sense[0] == 'Max' else c * x
        return lp, x, mip.A, mip.b, mip.sense[1], mip.integerIndices
    elif file_name is not None:
        lp = CyClpSimplex()
        m = lp.extractCyLPModel(file_name)
        x = m.getVarByName('x')
        integerIndices = [i for (i, j) in enumerate(lp.integerInformation) if j == True]
        infinity = lp.getCoinInfinity()
        sense = None
        for i in range(lp.nRows):
            if lp.constraintsLower[i] > -infinity:
                if sense == '<=':
                    print "Function does not support mixed constraint..."
                    break
                else: 
                    sense = '>='
                    b = lp.constraintsLower
            if lp.constraintsUpper[i] < infinity: 
                if sense == '>=':
                    print "Function does not support mixed constraint..."
                    break
                else: 
                    sense = '<='
                    b = lp.constraintsUpper
        return lp, x, lp.coefMatrix, b, sense, integerIndices
    else:
        print "No file or module name specified..."
        return None, None, None, None, None, None

if __name__ == '__main__':
    
    generate_bound_splits = True
    generate_max_viol_splits = False
    generate_GMI = False
    debug_print = True
    epsilon = 0.01
    maxiter = 1
    max_cuts = 10
    display = True
    if not DISPLAY_ENABLED:
        display = False
    
    lp, x, A, b, sense, integerIndices = read_instance(module_name = 'MIP6')
    #lp, x, A, b, sense, integerIndices = read_instance(file_name = 'p0033.mps')
    infinity = lp.getCoinInfinity()
    lp.logLevel = 0
    
    if display:
        disp_relaxation(A, b)
    
    for i in xrange(maxiter):
        print 'Iteration ', i
        lp.primal(startFinishOptions = 'x')
        print 'Current bound:', lp.objectiveValue
        bv = lp.basicVariables
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
        if isInt(sol[integerIndices]):
            print 'Integer solution found!'
            break
        cuts = []
        if generate_bound_splits:
            cuts += boundOptimalSplitCuts(lp, integerIndices, sense, max_coeff=10)
        if generate_max_viol_splits:
            cuts += maxViolationSplitCuts(lp, integerIndices, sense, max_coeff=10)
        if generate_GMI:
            cuts += gomoryCut(lp, integerIndices, sense)
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






