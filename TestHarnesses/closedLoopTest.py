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
import argparse
import random
import pickle
import traceback
import copy

PD = VCLC.PDControl()
PI = VCLC.PIControl()
PID = VCLC.PIDControl()

dT = 0.01
kp = 1.0
ki = 1.0
kd = 1.0
trim = 1.0
lowLimit = 0.0
highLimit = 10.0

command = 10.0
current = 1.0
derivative = 2.0

PD.setPDGains(kp,kd,trim,lowLimit,highLimit)
u = PD.Update(command,current,derivative)
print("PD u:" + str(u))

PI.setPIGains(dT,kp,ki,trim,lowLimit,highLimit)
u = PI.Update(command, current)
print("PI u:" + str(u))
PI.resetIntegrator()

PID.setPIDGains(dT,kp,kd,ki,trim,lowLimit,highLimit)
print("PID highLimit:" + str(PID.highLimit))
print("PID u:" + str(u))

controls = Controls.controlGains()

controls.kp_roll = 1.0
controls.kd_roll = 1.0
controls.ki_roll = 1.0

controls.kp_sideslip = 1.0
controls.ki_sideslip = 1.0

controls.kp_course = 1.0
controls.ki_course = 1.0

controls.kp_pitch = 1.0
controls.kd_pitch = 1.0
controls.kp_altitude = 1.0
controls.ki_altitude = 1.0
controls.kp_SpeedfromThrottle = 1.0
controls.kiSpeedfromThrottle = 1.0
controls.kp_SpeedfromElevator = 1.0
controls.ki_SpeedfromElevator = 1.0


controlLoop = VCLC.VehicleClosedLoopControl()
controlLoop.setControlGains(controls)
controlLoop.Update()
controls = controlLoop.getControlGains()
print(controls)
print(controlLoop.getVehicleAerodynamicsModel())
print(controlLoop.getVehicleControlSurface())
print(controlLoop.getVehicleState())
controlLoop.reset()
controlLoop.setTrimInputs(Inputs.controlInputs(0.7,0.3,0.7,0.1))
print(controlLoop.trimInputs.Throttle)
print(controlLoop.trimInputs.Aileron)
print(controlLoop.trimInputs.Elevator)
print(controlLoop.trimInputs.Rudder)

#controlLoop.setVehicleState()

