numVars = 2
points = [[2.5, 4.5], [6.5, 0.5], [0.5, 1],
          [7, 5.7], [7.7, 5], [2, 0.25]]
#points = [[0, 0], [2.5, 4.5], [7.75, 5.75], [6.5, 0.5]]
#points = [[0, 0], [2.5, 4.5], [7.5, 5.5], [6.5, 0.5]]
rays = []
c = [0, 1]
sense = ('Max', '<=')
integerIndices = [0, 1]

if __name__ == '__main__':

    try:
        from coinor.cuppy.cuttingPlanes import solve, gomoryCut
        from coinor.cuppy.milpInstance import MILPInstance
    except ImportError:
        from src.cuppy.cuttingPlanes import solve, gomoryCut
        from src.cuppy.milpInstance import MILPInstance

    m = MILPInstance(c = c, points = points, rays = rays, 
                     sense = sense, integerIndices = integerIndices,
                     numVars = numVars)
    
    solve(m, whichCuts = [(gomoryCut, {})], display = True, debug_print = True)
