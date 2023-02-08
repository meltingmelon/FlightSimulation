import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import math

import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Utilities.Rotations as Rotations
import ece163.Modeling.VehicleDynamicsModel as VDM
import ece163.Containers.States as States
import ece163.Modeling.VehicleAerodynamicsModel as VAM
from ece163.Containers import Inputs

print("Does this shit work?\n")
aero = VAM.VehicleAerodynamicsModel()

aero.vehicleDynamics.u = 10
aero.vehicleDynamics.v =11
aero.vehicleDynamics.w = 20
aero.wind.Wu = 2
aero.wind.Wv = 2
aero.wind.Ww = 2
aero.vehicleDynamics.Va = 10
aero.wind.Wn = 2
aero.wind.We = 1
aero.wind.Wd = 10

aero.CalculateAirspeed(aero.vehicleDynamics.state, aero.wind)
print("\n Print out the easy functions")
print(aero.getVehicleState())
print(aero.getWindState())

gravity = aero.gravityForces(aero.vehicleDynamics.state)
print("Gravity:" + str(gravity) +"\n")

aeroForces = aero.aeroForces(aero.vehicleDynamics.state)
print("AeroForces: " + str(aeroForces) +"\n")

controlInputs = Inputs.controlInputs()
propForces = aero.CalculatePropForces(aero.vehicleDynamics.state.Va, controlInputs.Throttle)
print ( "Prop Forces: " + str(propForces) + "\n")

control = aero.controlForces(controlInputs)
print("Control Forces: " + str(control) + "\n")

totalForces = aero.updateForces(aero.vehicleDynamics.state,  aero.wind, controlInputs)
print("Total Forces: " + str(totalForces) + "\n")


aero.Update(controlInputs)
print("Updated Aero:" + str(aero.getVehicleState()) + "\n")

aero.setWindModel()
