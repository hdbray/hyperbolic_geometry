import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import time
start_time = time.time()

### generators of PSL(2,Z)

A=[[1,1],[0,1]]
B=[[0,1],[1,0]]
Ainv=np.linalg.inv(A)
Binv=np.linalg.inv(B)

### function for matrix product 

def mp(list_of_matrices):
    X = np.identity(len(list_of_matrices[0]))
    for matrix in list_of_matrices:
        X = np.dot(X,matrix)
    return X.tolist()

### action of PSL(2,R) on upper half plane by mobius transformations
### this function is more general than needed for this picture but maybe it
### will be useful later?

def mobius_transformation(M,v):
    # M is a matrix which represents a mobius transformation
    # v is a vector in R^2 which will be our model for the complex plane
    # via the bijection [x,y] <-> x+iy
    a=M[0][0]
    b=M[0][1]
    c=M[1][0]
    d=M[1][1]
    x=v[0]
    y=v[1]
    numerator_real_part = (a*x+b)*(c*x+d)+a*c*y*y
    numerator_imaginary_part = (a*d-b*c-a*c*x)*y
    denomenator = (c*x+d)**2+(c*y)**2
    return [numerator_real_part/denomenator,numerator_imaginary_part/denomenator]

# need to upgrade this function so that it can also act on vertical
# geodesics

def mobius_of_geodesic(M,l,height='3'):
    # M is a matrix which represents a mobius transformation
    # l is an array two distinct vectors which represent end points of
    # geodesic in the upper half plane. If the vectors have the same
    # horizontal coordinate, that means the geodesic is a vertical geodesic 
#    a=M[0][0]
#    b=M[0][1]
#    c=M[1][0]
#    d=M[1][1]
    x_minus=l[0]
    x_plus=l[1]
    return [mobius_transformation(M,x_minus),mobius_transformation(M,x_plus)]

### generate the group

def group_elts(n,generators=[A,B,Ainv,Binv]):
    list_these_matrices=[ [] for i in range(n+1)]
    list_these_matrices[0]=[(np.identity(len(generators[0]))).tolist()]
    with_duplicate_matrices=[ [] for i in range(n+1)]
    with_duplicate_matrices[0]=[(np.identity(len(generators[0]))).tolist()]
    flat_list = []
    for j in range(n):
        for i in list_these_matrices[j]:
#        range(len(list_these_matrices[j])):
            for M in generators:
                with_duplicate_matrices[j+1].append(mp([M,i]))
                with_duplicate_matrices[j+1].append(mp([i,M]))
        for k in with_duplicate_matrices[j+1]:
            if k not in list_these_matrices[j+1]:
                list_these_matrices[j+1].append(k)
        flat_list.extend(list_these_matrices[j+1])
    return flat_list

#### draw some lines
#
#def draw_lines(list_of_vectors, this_color='black'):
#    for i in range(len(list_of_vectors)-1):
#        plt.plot([list_of_vectors[i][0],list_of_vectors[i+1][0]],[list_of_vectors[i][1],list_of_vectors[i+1][1]], color=this_color)
#
#
#draw_lines([[-3,0],[3,0]])
#for i in range(-2,3):
#    draw_lines([[i,0],[i,3]])

#n=3
#for M in group_elts(n):

patches.Arc((0,0), 1, 1, 0.0, 0.0, 359.9)
#matplotlib.patches.Arc(xy, width, height, angle=0.0, theta1=0.0, theta2=360.0,

axes = plt.gca()
axes.set_xlim([-2.2,2.2])
axes.set_ylim([-.2,2])
plt.show()
