numVars = 2
numCons = 5
A = [#[4, 1],
     [1, 4],
#     [1, -1],
     [-1, 0],
     [0, -1]] 

b = [#28,
     27,
#     1,
     0,
     0]

c = [2, 5]
obj_val = 2 
sense = ('Max', '<=')
integerIndices = [0, 1]
points = []
rays = []
rhs = None
cuts = None
