import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math
import ece163.Utilities.Rotations as Rotations
import ece163.Utilities.MatrixMath as MatrixMath

# Testing Limits of DCM2euler
DCM = [[0.5, 2, 2], [0.5, .1, .2], [0.2,1,3]]
yaw, pitch, roll = Rotations.dcm2Euler(DCM)
print(str(yaw) + "\n")
print(str(pitch) + "\n")
print(str(roll) + "\n")