numVars = 2
numCons = 6
A = [[ -8,     30],
     [ -2,    -4],
     [-14,     8],
     [  2,   -36],
     [  30, -8],
     [  10,     10],
#     [  1,    -8],
     ]
b = [115,
     -5,
      1,
     -5,
     191,
     127,
#     2
    ]
sense = ('Max', '<=')
integerIndices = [0, 1]
points = None
#Variations
#points = [[2.5, 4.5], [6.5, 0.5], [0.5, 1],
#          [7, 5.7], [7.7, 5], [2, 0.25]]
#points = [[0, 0], [2.5, 4.5], [7.75, 5.75], [6.5, 0.5]]
#points = [[0, 0], [2.5, 4.5], [7.5, 5.5], [6.5, 0.5]]
rays = []

#This is the integer hull
#A = [[ 0, -1],
# [-1, 1],
# [-2,  1],
# [-1,  2],
# [ 2, -1],
# [ 0,  1],
# [ 1, 0],
# ]

#b = [ -1,   
#   1,  
#  -1,   
#   5,  
#  11,   
#   5,   
#   7]

c = [-1, 1]
#c = [1, -1]
obj_val = 2 

#Rank 1 inequalities
cuts = [
     [-1.26, 2.72], 
#    pypolyhedron doesn't like the full-recision version 
#     [-1.25842697, 2.71910112], 
     [-2.71910112,  2.69662921],
     [ 2.06015038, -1.08270677], 
     [1, 0],
     [0, 1],
     [-1, 0],
     [0, -1],
     ]
rhs = [
     8.08,
#    pypolyhedron doesn't like the full-recision version 
#     8.08988764045,
     4.33707865169,
#     -0.563909774436,
     11.8496240602,
     7,
     5,
     -1,
     -1,
    ]

# This is to simulate branching
obj_val1 = 1.63
A1 = A + [[1, 0]]
b1 = b + [2]

A2 = A + [[-1, 0]]
b2 = b + [-3]

obj_val2 = 1
A3 = A2 + [[0, -1]]
b3 = b2 + [-5]

A4 = A2 + [[0, 1]]
b4 = b2 + [4]
