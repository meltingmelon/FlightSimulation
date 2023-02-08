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

import argparse
import random
import pickle
import traceback
import copy

linearizeModel = Linearized.transferFunctions()
tuning = Controls.controlTuning()

tuning.Wn_roll = 10
tuning.Zeta_roll = 1.0
tuning.Wn_course = 1.0	# Wn_roll should be 5-10x larger
tuning.Zeta_course = 1.0
tuning.Wn_sideslip = 1.0
tuning.Zeta_sideslip = 1.0
#tuning knows for longitudinal control
tuning.Wn_pitch = 10.0
tuning.Zeta_pitch = 10.0
tuning.Wn_altitude = 6.0	# Wn_pitch should be 5-10x larger
tuning.Zeta_altitude = 1.0
tuning.Wn_SpeedfromThrottle = 1
tuning.Zeta_SpeedfromThrottle = 5.0
tuning.Wn_SpeedfromElevator = 3.0
tuning.Zeta_SpeedfromElevator = 0.5

print(tuning)

linearizeModel.Va_trim = 2.0
linearizeModel.alpha_trim = 0.4
linearizeModel.beta_trim = 0.5
linearizeModel.gamma_trim = 0.4
linearizeModel.theta_trim = 0.0
linearizeModel.phi_trim = 0.4
linearizeModel.a_phi1 = 0.5
linearizeModel.a_phi2 = 1.0
linearizeModel.a_beta1 = 1.0
linearizeModel.a_beta2 = 1.0
linearizeModel.a_theta1 = 2.0
linearizeModel.a_theta2 = 0.5
linearizeModel.a_theta3 = 1.0
linearizeModel.a_V1 = 10
linearizeModel.a_V2 = 1.0
linearizeModel.a_V3 = 0.5

print("Does this shit work")
control = Controls.controlGains()
control = VCG.computeGain(tuning, linearizeModel)
print(control)
print("computeGains works\n")

tuning = VCG.computeTuningParameters(control, linearizeModel)
print(tuning)
print("computeTuningParameters works")



