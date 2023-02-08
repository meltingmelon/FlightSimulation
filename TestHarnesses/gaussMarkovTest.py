import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math

import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Utilities.Rotations as Rotations
import ece163.Modeling.VehicleAerodynamicsModel as VAM
import ece163.Containers.Inputs as Inputs
import ece163.Containers.States as States
import ece163.Containers.Linearized as Linearized
import ece163.Containers.Controls as Controls
import ece163.Modeling.WindModel as WindModel
import ece163.Controls.VehiclePerturbationModels as VPM
import ece163.Controls.VehicleTrim as VT
import ece163.Controls.VehicleControlGains as VCG
import ece163.Controls.VehicleClosedLoopControl as VCLC
import ece163.Sensors.SensorsModel as SM
import argparse
import random
import pickle
import traceback
import copy


gm = SM.GaussMarkov(0.01, 1e6, 1.0)
gmXYZ = SM.GaussMarkovXYZ(0.01, 1e6,1.0)

print("v default: ", gm.vnoise)
gm.update()
print("v update: ", gm.vnoise)
gm.reset()
print("v reset: ", gm.vnoise)

print("X default: ", gmXYZ.vNoiseX.vnoise)
print("Y default: ", gmXYZ.vNoiseY.vnoise)
print("Z default: ", gmXYZ.vNoiseZ.vnoise)

x,y,z = gmXYZ.update()

print("X updated: ", x)
print("Y updated: ", y)
print("Z updated: ", z)

gmXYZ.reset()

print("X reset: ", gmXYZ.vNoiseX.vnoise)
print("Y reset: ", gmXYZ.vNoiseY.vnoise)
print("Z reset: ", gmXYZ.vNoiseZ.vnoise)


