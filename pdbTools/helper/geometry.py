# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
geometry.py

A set of geometry functions for manipulating pdb files.
"""
__author__ = "Michael J. Harms"
__date__ = "080201"

from math import sqrt, cos, sin, acos, pi

def dist(c1,c2=[0.,0.,0.]):
    """
    Calculate the distance between two coordinates in 3d space.  If no c2 is
    specified it returns the length of the vector.
    """
        
    return sqrt(sum([(c1[i] - c2[i])**2 for i in range(3)]))

def dist_sq(c1,c2=[0.,0.,0.]):
    """
    Calculate the squared distance between two coordinates in 3d space.  If no
    c2 is specified it returns the squared length of the vector.
    """
        
    return sum([(c1[i] - c2[i])**2 for i in range(3)])

def calcDistances(coord):
    """
    Calculate all distances in coord.
    """
    
    num_points = len(coord)
    d = [[0. for j in range(num_points)] for i in range(num_points)]
    for i in xrange(num_points):
        d[i][i] = 0.
        for j in xrange(i+1,num_points):
            d[i][j] = dist(coord[i],coord[j])
            d[j][i] = d[i][j]

    return d

def crossProduct(u,v):
    """
    Calculates the cross product of two 3d vectors (as 1-d arrays).
    """

    prod = [0.,0.,0.]
    prod[0] = u[1]*v[2]-u[2]*v[1]
    prod[1] = u[2]*v[0]-u[0]*v[2]
    prod[2] = u[0]*v[1]-u[1]*v[0]
    
    return prod

def dotProduct(u,v):
    """
    Calculates the dot product between two vectors.
    """

    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def findAngle(u,v):
    """
    Calculates the angle (degrees) between two vectors (as 1-d arrays) using
    dot product.
    """
    
    mag_u = sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    mag_v = sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    
    return 180/pi * acos(dotProduct(u,v)/(mag_u*mag_v))

def genRotMatrix(axis,theta):
    """
    Generate a rotation matrix for rotation of theta about axis.
    """

    matrix = [[0. for j in range(3)] for i in range(3)]

    axis_length = sqrt((axis[0]**2 + axis[1]**2 + axis[2]**2))
    xNorm = axis[0]/axis_length
    yNorm = axis[1]/axis_length
    zNorm = axis[2]/axis_length

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    one_costheta = 1.0 - cos_theta

    matrix[0][0] = cos_theta                + xNorm*xNorm*one_costheta
    matrix[0][1] = xNorm*yNorm*one_costheta - zNorm*sin_theta
    matrix[0][2] = xNorm*zNorm*one_costheta + yNorm*sin_theta
    matrix[1][0] = xNorm*yNorm*one_costheta + zNorm*sin_theta
    matrix[1][1] = cos_theta                + yNorm*yNorm*one_costheta
    matrix[1][2] = yNorm*zNorm*one_costheta - xNorm*sin_theta
    matrix[2][0] = xNorm*zNorm*one_costheta - yNorm*sin_theta
    matrix[2][1] = yNorm*zNorm*one_costheta + xNorm*sin_theta
    matrix[2][2] = cos_theta                + zNorm*zNorm*one_costheta

    return matrix


def arbRot(vector,axis,theta):
    """
    Rotates a vector about an arbitrary axis by an arbitrary amount.
    """

    matrix = genRotMatrix(axis,theta)   
    
    return [sum([vector[j]*matrix[j][i] for j in range(3)]) for i in range(3)]


def arbRotCoord(coord,axis,theta):
    """
    Rotate all vectors in coord about an arbitray axis by theta.
    """

    matrix = genRotMatrix(axis,theta)
    
    return [[sum([c[j]*matrix[j][i] for j in range(3)]) for i in range(3)]
            for c in coord]
    


def calcGlyCbeta(Ncoord,CAcoord,COcoord):
    """
    Generates a beta carbon for a glycine using the coordinates for the amide N,
    alpha C, and carboxyl C.
    """

    CA_CO_vector = []; CA_N_vector = []
    for i in range(3):
        CA_CO_vector.append(COcoord[i] - CAcoord[i])
        CA_N_vector.append(Ncoord[i] - CAcoord[i])

    rotation_amount = 240*(pi/180.)

    rotated = arbRot(CA_CO_vector, CA_N_vector, rotation_amount)

    CBeta = []
    for i in range(3):
        CBeta.append(rotated[i] + CAcoord[i])

    return CBeta   

 
def calcHXT(C_coord,O_coord,OXT_coord):
    """
    Calculates the location of HXT using the location of C, O, and OXT.
    (C-terminal hydrogen).
    """
        
    vector_1,vector_2,vector_3 = [0.,0.,0.], [0.,0.,0.], [0.,0.,0.]
    for i in range(3):
        vector_1[i] = O_coord[i] - C_coord[i]
        vector_2[i] = OXT_coord[i] - C_coord[i]
        vector_3[i] = vector_1[i] + vector_2[i]

    dv3 = sqrt(sum([vector_3[i]*vector_3[i] for i in range(3)]))

    HXT_coord = [0.,0.,0.]
    for i in range(3):
        vector_3[i] = vector_3[i]/dv3
        HXT_coord[i] = OXT_coord[i] + vector_3[i]

    return HXT_coord

def calcHG(CB_coord,SG_coord):
    """
    Calculates the location of HG using the location of CB and SG.
    (Hydrogen on free cysteines).
    """

    vector_1 = [0.,0.,0.]
    vector_1 = [SG_coord[i] - CB_coord[i] for i in range(3)]

    dv3 = sqrt(sum([vector_1[i]*vector_1[i] for i in range(3)]))

    HG_coord = [0.,0.,0.]
    for i in range(3):
        vector_1[i] = vector_1[i]/dv3
        HG_coord[i] = SG_coord[i] + 1.08*vector_1[i]

    return HG_coord

def calcDihedrals(prevCO,currN,currCA,currCO,nextN,cutoff=6.5):
    """
    Calculates phi and psi angles for an individual residue.
    """
    
    # Set CA coordinates to origin
    A = [prevCO[i] - currCA[i] for i in range(3)]
    B = [currN[i] - currCA[i] for i in range(3)]
    C = [currCO[i] - currCA[i] for i in range(3)]
    D = [nextN[i] - currCA[i] for i in range(3)]
 
    # Make sure the atoms are close enough 
    #if max([dist_sq(x) for x in [A,B,C,D]]) > cutoff:
    #    err = "Atoms too far apart to be bonded!"
    #    raise ValueError(err)
 
    # Calculate necessary cross products (define vectors normal to planes)
    V1 = crossProduct(A,B)
    V2 = crossProduct(C,B)
    V3 = crossProduct(C,D)

    # Determine scalar angle between normal vectors
    phi = findAngle(V1,V2)
    if dotProduct(A,V2) > 0: phi = -phi
    
    psi = findAngle(V2,V3)
    if dotProduct(D,V2) < 0: psi = -psi
    
    return phi, psi


def calcHN(prevCO,prevO,currN,currCA):
    """
    Calculate the position of the amide hydrogen.
    """

    CO_bond = [prevO[i] - prevCO[i] for i in range(3)]
    CN_bond = [currN[i] - prevCO[i] for i in range(3)]

    return [prevCO[i] + CO_bond[i] + CN_bond[i] for i in range(3)]

