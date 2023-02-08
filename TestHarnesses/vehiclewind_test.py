import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math
import ece163.Utilities.Rotations as Rotations
import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Modeling.WindModel as WindModel
import ece163.Containers.Inputs as Inputs
import ece163.Containers.States as States
import ece163.Modeling.VehicleDynamicsModel as VDM
import ece163.Constants.VehiclePhysicalConstants as VPC

print("Does this shit work?")

windModel = WindModel.WindModel()
windModel.CreateDydenTransferFns(windModel.dT, windModel.Va, windModel.drydenParamters)


print("PhiU: " + str(windModel.PhiU) + "\n")
print("GammaU: "+ str(windModel.GammaU) + "\n")
print("Hu:" + str(windModel.Hu) + "\n")

print("PhiV: " + str(windModel.PhiV) + "\n")
print("GammaV: " + str(windModel.GammaV) + "\n")
print("Hv: " + str(windModel.Hv) + "\n")

print("PhiW: " + str(windModel.PhiW) + "\n")
print("GammaW: " + str(windModel.GammaW) + "\n")
print("Hw: " + str(windModel.Hw) + "\n")

print(windModel.getWind())
windModel.reset()
wind = States.windState(1,2,2,1,2,1)
windModel.setWind(wind)
print(windModel.getWind())
windModel.reset()
print(windModel.getWind())

windModel.Update(0.2, 0.1, .6)
print("First Wind update")
print("WindU:" + str(windModel.Wind.Wu) + "\n")
print("WindV:" + str(windModel.Wind.Wv) + "\n")
print("WindW:" + str(windModel.Wind.Ww) + "\n")

windModel.Update(0.1, 0.2, .5)
print("Second Wind update")
print("WindU:" + str(windModel.Wind.Wu) + "\n")
print("WindV:" + str(windModel.Wind.Wv) + "\n")
print("WindW:" + str(windModel.Wind.Ww) + "\n")
