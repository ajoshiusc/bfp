# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 15:28:12 2016

@author: ajoshi
"""

from mayavi import mlab
import numpy as np
# from scipy.special import sph_harm
import math


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


# Create a sphere
r = 0.3
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:501j, 0:2 * pi:501j]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
mlab.clf()
# Represent spherical harmonics on the surface of the sphere
mlab.mesh(x, y, z, color=(0.5, 0.5, 0.5), opacity=1)

# phi, theta = np.mgrid[-.2:.2:10j, 0:2*pi:10j]
phi = 0.1 * np.random.randn(100)
theta = np.random.randn(100) * 2.0 * pi;
theta = theta + pi / 6.0
x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

R = rotation_matrix([2, 3, 1], pi / 4)
v = np.dot(R, np.array([x, y, z]))
mlab.points3d(v[0,], v[1,], v[2,], scale_factor=0.01, color=(1, 0, 0))

phi = 0.1 * np.random.randn(100)
theta = np.random.randn(100) * 2.0 * pi;
theta = theta + pi / 6.0
x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

R = rotation_matrix([-2, 3, 1], pi / 4)
v = np.dot(R, np.array([x, y, z]))
mlab.points3d(v[0,], v[1,], v[2,], scale_factor=0.01, color=(0, 1, 0))

phi = 0.1 * np.random.randn(100)
theta = np.random.randn(100) * 2.0 * pi;
theta = theta + pi / 6.0
x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

R = rotation_matrix([2, 3, 1], -pi / 4)
v = np.dot(R, np.array([x, y, z]))
mlab.points3d(v[0,], v[1,], v[2,], scale_factor=0.01, color=(0, 0, 1))
mlab.draw()
mlab.show(stop=True)
mlab.close()