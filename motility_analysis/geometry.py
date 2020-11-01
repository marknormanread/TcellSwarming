"""
Created on 28/10/2014

@author: Mark N. Read, 2017
"""
import math
import numpy


def radians_to_degrees(rad):
    return rad * 180.0 / numpy.pi


def degrees_to_radians(deg):
    return deg * numpy.pi / 180.0


def vector_length(v):
    """  Calculates the length of the vector with the supplied x y and z components.  """
    return math.sqrt( vector_length2(v) )


def vector_length2(v):
    """ Return the length squared of the supplied vector. """
    return sum(n * n for n in v)


def distance_between_points(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2
    return math.sqrt((dx * dx) + (dy * dy) + (dz * dz))


def dot_product(v, w):
    """
    Calculates the dot product between two vectors. The first vector has components x1, y1 and z1, the second x2, y2,
    and z2.
    """
    x1, y1, z1 = v
    x2, y2, z2 = w
    return (x1 * x2) + (y1 * y2) + (z1 * z2)


def cross_product(u, v):
    """ u and v should be tuples of (x, y, z). """
    u1, u2, u3 = u
    v1, v2, v3 = v
    s1 = u2 * v3 - u3 * v2
    s2 = u3 * v1 - u1 * v3
    s3 = u1 * v2 - u2 * v1
    return s1, s2, s3


class InvalidAngleException(Exception):
    """ Returned when a method is asked to compute and angle with an impossible solution. """
    pass


def angle_between_vectors(v, w):
    """
    A.B = |A||B|cos(alpha). Rearrange to give (A.B)/(|A||B|) = cos(alpha). Angle given in degrees.

    Raises InvalidAngleException if the two vectors have zero length.
    """
    dot = dot_product(v, w)
    lenV = vector_length(v)
    lenW = vector_length(w)
    if (lenV * lenW) == 0.0: 			# avoid div by zero.
        raise InvalidAngleException

    tmp = dot / (lenV * lenW)
    # necessary because of precision errors. tmp occasionally <-1 or >1, likely because values are bing generated
    #  in another program, and written with limited precision to files that are used for further calculation
    # here.
    tmp = max(-1.0, min(tmp,1))
    radians = math.acos(tmp)
    degrees = radians * (180 / math.pi)
    return degrees


def vector_projection(v, w):
    """
    Projects vector V onto vector W.

    projection = b (a.b) / (b.b)      - where (x.y) signifies dot product of two vectors.
    """
    w1, w2, w3 = w
    dotvw = dot_product(v, w)
    dotww = dot_product(w, w)
    length = dotvw / dotww
    return (length*w1, length*w2, length*w3)


#######################################################################################################################
### THE FOLLOWING IS QUATERNION BASED MATHS, TAKEN FROM
### http://stackoverflow.com/questions/4870393/rotating-coordinate-system-via-a-quaternion
###
### accessed on 28/10/2014. It has been modified in places.
#######################################################################################################################

def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        if mag2 == 0.0:
            raise Exception('div by zero, cant normalise a vector of length zero')
        mag = math.sqrt(mag2)
        v = tuple(n / mag for n in v)
        #print vs
    return v


def q_mult(q1, q2):
    """
    This is used to multiply rotations together. Rotations are represented as unit vectors, and the result of this
    is also a unit vector.
    """
    w1, x1, y1, z1 = q1			# tuple expansion.
    w2, x2, y2, z2 = q2
    w = w1 * w2  -  x1 * x2  -  y1 * y2  -  z1 * z2
    x = w1 * x2  +  x1 * w2  +  y1 * z2  -  z1 * y2
    y = w1 * y2  +  y1 * w2  +  z1 * x2  -  x1 * z2
    z = w1 * z2  +  z1 * w2  +  x1 * y2  -  y1 * x2
    return w, x, y, z


def q_conjugate(q):
    q = normalize(q)
    w, x, y, z = q
    return (w, -x, -y, -z)


def qv_mult(q1, v1):
    """
    Multiply a quaternion by a vector. q1 is a quaternion (len(q1) = 4), v1 is a vector (len(v1) = 3; ie, no w value).
    """
    q2 = (0.0,) + v1			# add the w value here, hence converting vector to quaternion.
    # perform q1 . v1 . q1-1 (conjugate)
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]


def qv_transform(q, v):
    """ Applies the transformation represented in q to vector v. """
    vecX, vecY, vecZ = v
    w, x, y, z = q
    new_x =   w*w*vecX + 2*y*w*vecZ - 2*z*w*vecY +   x*x*vecX + 2*y*x*vecY + 2*z*x*vecZ -   z*z*vecX - y*y*vecX
    new_y = 2*x*y*vecX +   y*y*vecY + 2*z*y*vecZ + 2*w*z*vecX -   z*z*vecY +   w*w*vecY - 2*x*w*vecZ - x*x*vecY
    new_z = 2*x*z*vecX + 2*y*z*vecY +   z*z*vecZ - 2*w*y*vecX -   y*y*vecZ + 2*w*x*vecY -   x*x*vecZ + w*w*vecZ
    return (new_x, new_y, new_z)


def axisangle_to_q(v, theta):
    v = normalize(v)
    x, y, z = v
    theta /= 2
    w = math.cos(theta)
    x = x * math.sin(theta)
    y = y * math.sin(theta)
    z = z * math.sin(theta)
    return w, x, y, z


def q_to_axisangle(q):
    w, v = q[0], q[1:]
    theta = math.acos(w) * 2.0
    return normalize(v), theta


########## USED FOR DEBUGGING PURPOSES #############
if __name__ == "__main__":
    x_axis_unit = (1, 0, 0)
    y_axis_unit = (0, 1, 0)
    z_axis_unit = (0, 0, 1)

    v = (1, 1, 1)
    w = (5,0,0)

    print(vector_projection(v, w))

