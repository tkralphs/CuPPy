numVars = 2
numCons = 3

# For display, adding bounds explicitly is better
A = [[-2, -3], 
     [-4, -2], 
     [-3, -4], 
     [1, 0], 
     [0, 1], 
     [-1, 0], 
     [0, -1]]

b = [-5, 
     -15, 
     -20, 
     9, 
     6, 
     0, 
     0]

#original data
'''
A = [[-2, -3], 
     [-4, -2], 
     [-3, -4],
     ]

b = [-5, 
     -15, 
     -20,
      ]

x_u = [9,
       6,
       ]

'''

integerIndices = [0, 1]

c = [20,
     15,
     ]

obj_val = 100

sense = ('Min', '<=')

cuts = [[-4.8, -4.4], 
        [-1, -1],
        [-2, -1],
        [-3, -2],
        [-2, -1],
        ]

rhs = [-26, 
       -6, 
       -8,
       -14,
       -10, 
       ]

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
