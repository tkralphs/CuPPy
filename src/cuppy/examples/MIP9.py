numVars = 2
numCons = 4
#points = [[0, 0], [3, 4], [8, 6], [6, 1]]
A = [[-4, 3],
     [1, -6],
     [5, -2],
     [-2, 5]]
b = [0, 0, 28, 14]
sense = ('Max', '<=')
integerIndices = [0, 1]
c = [1, 1]

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
