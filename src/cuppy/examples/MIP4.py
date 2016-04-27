numVars = 2
numCons = 4
A = [[1, -1],
     [-1, 3],
     [-7, 3],
     [3, -7],
     ]
b = [0, 0, -18, -8]
c = [-1, -1]
sense = ('Max', '<=')
integerIndices = [0, 1]

if __name__ == '__main__':

    try:
        from coinor.cuppy.cuttingPlanes import solve, gomoryCut
        from coinor.cuppy.milpInstance import MILPInstance
    except ImportError:
        from src.cuppy.cuttingPlanes import solve, gomoryCut
        from src.cuppy.milpInstance import MILPInstance

    m = MILPInstance(A = A, b = b, c = c, 
                     sense = sense, integerIndices = integerIndices,
                     numVars = numVars)
    
    solve(m, whichCuts = [(gomoryCut, {})], display = True, debug_print = True)
