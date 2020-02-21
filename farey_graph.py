import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import time
start_time = time.time()

### generators of PSL(2,Z)

A=[[1,1],[0,1]]
B=[[0,1],[1,0]]
Ainv=np.linalg.inv(A)
# note that B is its own inverse so we do not need to introduce it 

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
    # output will be a vector
    a=M[0][0]
    b=M[0][1]
    c=M[1][0]
    d=M[1][1]
    x=v[0]
    y=v[1]
    numerator_real_part = (a*x+b)*(c*x+d)+a*c*y*y
    numerator_imaginary_part = (a*d-b*c-a*c*x)*y
    denominator = (c*x+d)**2+(c*y)**2
#    if denominator!=0:
    return [numerator_real_part/denominator,numerator_imaginary_part/denominator]
#    else:

# need to upgrade this function so that it can also act on vertical
# geodesics

def mobius_of_geodesic(M,r,t,height=4):
    # M is a matrix which represents a mobius transformation
    # a,b are distinct vectors with non-negative vertical coordinate. 
    # the function should calculate the geodesic from a to b and return the
    # pair of endpoints at infinity of the geodesic - except if the
    # geodesic is a vertical line, then one of the endpoints will have
    # vertical coordinate equal to the height (this will represent
    # infinity) 


    # for the moment I am going to assume that if a,b have different
    # horizontal coordinates, then a,b are both on
    # the real axis. It's much simpler. In the long run, would be
    # nice to have a function which takes any two points and
    # calculates their geodesic orthocircle
    r_hor=r[0]
    r_vert=r[1]
    t_hor=t[0]
    t_vert=t[1]
    a=M[0][0]
    b=M[0][1]
    c=M[1][0]
    d=M[1][1]
    if r_hor==t_hor:
        if c==0:
            endpoint1=mobius_transformation(M,[r_hor,0])
            endpoint2=mobius_transformation(M,[r_hor,height])
        else: 
            if c*r_hor+d==0:
                endpoint1=[a/c,height]#mobius_transformation(M,[r_hor,0])
                endpoint2=[a/c,0.]
            else:
                endpoint1=mobius_transformation(M,[r_hor,0])
                endpoint2=[a/c,0.]
    else:
        endpoint1=mobius_transformation(M,r)
        endpoint2=mobius_transformation(M,t)
    return [endpoint1,endpoint2]


#print(mobius_of_geodesic(mp([A,B]),[1,3],[1,0]))

### generate the group

def group_elts(n,generators=[A,B,Ainv]):
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
def draw_lines(list_of_vectors, default_color='black',default_linewidth=1):
    for i in range(len(list_of_vectors)-1):
        plt.plot([list_of_vectors[i][0],list_of_vectors[i+1][0]],[list_of_vectors[i][1],list_of_vectors[i+1][1]], color=default_color, linewidth=default_linewidth)
#
#
#draw_lines([[-3,0],[3,0]])
#for i in range(-2,3):
#    draw_lines([[i,0],[i,3]])

#n=3
#for M in group_elts(n):

##### draw a geodesic orthocircle
def orthocircle(a,b,default_color='black',default_linewidth=1):
    # a and b are real numbers
    center=(a+b)/2
    radius=abs(a-b)
    return patches.Arc((center,0), radius, radius, 0.0, 0.0, 180, color=default_color, linewidth=default_linewidth)

#geodesic= patches.Arc((0,0), 1, 1, 0.0, 0.0, 180)
#patches.Arc((xcenter,ycenter), width, height, angle=0.0, theta1=0.0,
#theta2=360.0)

def draw_geodesic(vectors, default_linewidth=1):
    v=vectors[0]
    w=vectors[1]
    v_hor=v[0]
    v_vert=abs(v[1])
    w_hor=w[0]
    w_vert=abs(w[1])
    if v_hor==w_hor:
#        draw_lines([v,w])
        draw_lines([[v_hor,v_vert],[w_hor,w_vert]],'black',default_linewidth)
    else:
        ax.add_patch(orthocircle(v_hor,w_hor,'black',default_linewidth))
    

fig = plt.figure()
ax = fig.add_subplot()

#ax.add_patch(orthocircle(-3,1))
#ax.add_patch(orthocircle(-1,1))

lw=.3

draw_lines([[-3,0],[3,0]])
draw_lines([[0,4],[0,0]],'black',lw)

n=15
for G in group_elts(n):
    draw_geodesic(mobius_of_geodesic(G, [0,4],[0,0]),lw)

#print(group_elts(n))

#draw_geodesic(mobius_of_geodesic(A, [0,4],[0,0]))
#draw_geodesic(mobius_of_geodesic(Ainv, [0,4],[0,0]))
#draw_geodesic(mobius_of_geodesic(B, [0,4],[0,0]))


#draw_geodesic([[1,4],[1,0]])
#draw_geodesic([[2,0],[1,0]])

axes = plt.gca()
axes.set_xlim([-2.2,2.2])
axes.set_ylim([-.2,2.2])

ax.set_aspect("equal")


plt.show()
