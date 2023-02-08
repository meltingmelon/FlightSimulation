import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math
from ece163.Modeling import VehicleAerodynamicsModel  as VAD
from ece163.Constants import VehiclePhysicalConstants as VPC
from ece163.Containers import States
from ece163.Containers import Inputs
from ece163.Containers import Linearized
from ece163.Utilities import MatrixMath
from ece163.Utilities import Rotations
from ece163.Controls import VehicleTrim
from ece163.Controls import VehiclePerturbationModels as VPM
from ece163.Modeling import WindModel


Va = 1
Throttle = 1
dThdThr = VPM.dThrust_dThrottle(Va, Throttle)
print("d Thrust dThrottle:" + str(dThdThr) + "\n")

dThedVa = VPM.dThrust_dVa(Va, Throttle)
print("d thrust dVa:" + str(dThedVa) + "\n")

vTrim = VehicleTrim.VehicleTrim()
Vastar = 25.0
Gammastar = math.radians(6.0)
Kappastar = -1.0 / 150.0

transfer = Linearized.transferFunctions()

check = vTrim.computeTrim(Vastar, Kappastar, Gammastar)
if check:
 print("Optimization successful")
else:
 print("Model converged outside of valid inputs, change parameters and try again")

wind = WindModel.WindModel()
vad = VAD.VehicleAerodynamicsModel()
vad.CalculateAirspeed(vad.vehicleDynamics.state,wind.Wind )


transfer = VPM.CreateTransferFunction(vTrim.VehicleTrimModel.vehicleDynamics.state, vTrim.ControlTrim)
print("Transfer Functions:  " + str(transfer) + "\n")
print("Va_trim: " + str(transfer.Va_trim) + "\n")
print("alpha_trim: " + str(transfer.alpha_trim)+ "\n")