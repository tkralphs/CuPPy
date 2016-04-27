numVars = 2
numCons = 4

A = [
     [1, 1],
     [4, -10],
     [-2, -2],
     [-6, -2],
     [-1, 4],
     [-7, 1],
     [0, -1],
     [1, -1],
     [4, 1],
     [0, 1],
     [-1, 5]
     ]
b = [
     8,
     -3,
     -9,
     -19,
     12,
     -13,
     -1,
     3,
     27,
     5,
     20
     ]

#points = [[0, 0], [3, 4], [8, 6], [6, 1]]

sense = ('Min', '<=')
integerIndices = [0, 1]
c = [1, 0]

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
