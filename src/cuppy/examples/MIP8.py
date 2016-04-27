numVars = 2
A = [[ 1, -2],
 [-0, -1],
 [-2, -1],
 [-1, -0],
 [-1, -1],
 [1, -3],
 [0,  1],
 [1, 0]]

b = [4,
     0,   
    -3,
    0,
    -2,
    3,
    4,
    8]

c = [-3, 5]
sense = ['Max', '<=']
integerIndices = [0, 1]

#This is for the separation code
x = [1, 1]
x_sep = [0, 1]

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
