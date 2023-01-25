import numpy as np

def calc_dihedral(a, b, c, d):
    ab = b - a 
    bc = c - b
    cd = d - c

    n1 = np.cross(ab, bc)
    n2 = np.cross(bc, cd)

    x = np.dot(n1, n2)
    y = np.dot(np.cross(bc, n1), n2)
    return np.sign(np.arctan2(y, x)) * np.arctan2(np.linalg.norm( np.cross(n1,n2) ) , np.dot(n1, n2))


def dist(a, b):
    a = np.array(a)
    b = np.array(b)
    return np.linalg.norm(a-b)


def calc_angle_b(a, b, c):
    if np.all(b == c):
        return np.pi / 2

    C = dist(a, b)
    A = -1*dist(b, c)
    B = dist(a, c)
    angle = np.arccos((np.dot(b-a,c-b)/(C*A)))
    if angle > np.pi:
        return np.pi - angle
    return angle

