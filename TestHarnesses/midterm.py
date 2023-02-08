import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math

import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Utilities.Rotations as Rotations
import ece163.Modeling.VehicleAerodynamicsModel as VAM
import ece163.Modeling.WindModel as WM
import ece163.Containers.States as States
import ece163.Containers.Inputs as Inputs
import ece163.Containers.Controls as Controls
import ece163.Controls.VehicleControlGains as CG
import ece163.Controls.VehicleClosedLoopControl as CLC
import argparse
import random
import pickle
import traceback
import copy
import numpy as np


#Euler angles
yaw = 0.6
pitch = 0.5
roll = 2.0

va = 2.0
vg = 3.0
wn = 0.5
we = 0.3
wd = 0.3

alpha = 0.2 # angle of attack
beta = 0.1 # sideslip angle,

gamma = 0.4 # flight path angle
chi = 0.3 # course angle
crab = chi - yaw # crab angle
gamma_airmass = pitch - alpha

u = 0.2
v = 0.5
w = 0.9

# Airspeed vector w.r.t. intertial frame
va_rot = [[math.cos(yaw) * math.cos(gamma_airmass)], [math.sin(yaw) * math.cos(gamma_airmass)], [-math.sin(gamma_airmass)]]
va_i = MatrixMath.matrixScalarMultiply(va, va_rot)

# DCM
DCM = [[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]

# Euler to DCM conversion
DCM_out = Rotations.euler2DCM(yaw, pitch, roll)
print("DCM: ", DCM_out)

# DCM to euler
yaw_out, pitch_out, roll_out = Rotations.dcm2Euler(DCM)
print("yaw_out: ", yaw_out)
print("pitch_out: ", pitch_out)
print("roll_out: ", roll_out)

#position vector
p = [[.5],[.3],[1.3]]

# Rotation matrices
# v -> v1 
R_v1 = [[math.cos(yaw), math.sin(yaw), 0.0],
        [-math.sin(yaw), math.cos(yaw), 0.0],
        [0.0, 0.0, 1.0]]

# v1 -> v2
R_v2 = [[math.cos(pitch), 0.0, -math.sin(pitch)],
        [0.0, 1.0, 0.0],
        [math.sin(pitch), 0.0, math.cos(pitch)]]

# v2 -> b, Use DCM_out for i -> b
R_b = [[1.0, 0.0, 0.0],
        [0.0, math.cos(roll), math.sin(pitch)],
        [0.0, -math.sin(pitch), math.cos(pitch)]]

# b -> s:
R_s = [[math.cos(alpha), 0.0, math.sin(alpha)],
        [0.0, 1.0, 0.0],
        [-math.sin(alpha), 0.0, math.cos(alpha)]]

# s -> w
R_w = [[math.cos(beta), math.sin(beta), 0.0],
        [-math.sin(beta), math.cos(beta), 0.0],
        [0.0, 0.0, 1.0]]

# b -> w
R_bw = MatrixMath.matrixMultiply(R_s,R_w)


# Some kind of rotation operation
rotated_position = MatrixMath.matrixMultiply(R_bw,p)
print("Rotated coords: ", rotated_position)

# Coordinate Turn
m = 1.0
g = 9.8
f_lift = (m * g) * (math.cos(gamma) / math.cos(roll))

#vg = va + vw
vg = 1.0
va = 1.0
vw = 1.0
chi_dot = (g * math.tan(roll) / (vg))

radius = (vg * math.cos(gamma) / chi_dot)
yaw_rate = (g / va) * math.tan(roll)

