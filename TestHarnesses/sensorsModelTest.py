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

sensors = SM.SensorsModel()
noisySensors = sensors.getSensorNoisy()
sensorsTrue = sensors.getSensorsTrue()
print("SensorsNoisy (before update): ", noisySensors)
print("SensorsTrue (before update): ", sensorsTrue)
print("SensorsBiases: ", sensors.sensorBiases)
print("SensorsSigmas: ", sensors.sensorSigmas)

print("Gyro Gauss Markov: ", sensors.gyroGM)
print("GPS Gauss Markove: ", sensors.gpsGM)

state = States.vehicleState(1,0.4,0.5,1,2,0.4,0.5,0.8,0.1,1,2,.7,)
sensors.VAM.vehicleDynamics.setVehicleState(state)

dot = States.vehicleState(1,1,1,1,1,1,1,1,1,1,1,1)
sensors.VAM.setVehicleDerivative(dot)

accel_x,accel_y,accel_z = sensors.updateAccelsTrue(sensors.VAM.vehicleDynamics.state, sensors.VAM.vehicleDynamics.dot)
print("accel x", accel_x)
print("accel_y", accel_y)
print("accel_z", accel_z)

gps_n, gps_e, gps_alt, gps_sog, gps_cog = sensors.updateGPSTrue(sensors.VAM.vehicleDynamics.state, sensors.VAM.vehicleDynamics.dot)
print("gps_n", gps_n)
print("gps_e", gps_e)
print("gps_alt", gps_alt)
print("gps_sog", gps_sog)
print("gps_cog", gps_cog)

gyro_x, gyro_y, gyro_z = sensors.updateGyrosTrue(sensors.VAM.vehicleDynamics.state)
print("gyro_x", gyro_x)
print("gyro_y", gyro_y)
print("gyro_Z", gyro_z)

mag_x,mag_y,mag_z = sensors.updateMagsTrue(sensors.VAM.vehicleDynamics.state)
print("mag_x", mag_x)
print("mag_y", mag_y)
print("mag_z", mag_z)

baro, pitot = sensors.updatePressureSensorsTrue(sensors.VAM.vehicleDynamics.state)
print("baro", baro)
print("pitot", pitot)

sensors.sensorsTrue = sensors.updateSensorsTrue(sensors.sensorsTrue,state,dot)
print("sensors True updated: ", sensors.sensorsTrue)

sensors.sensorsNoisy = sensors.updateSensorsNoisy(sensors.sensorsTrue,sensors.sensorsNoisy, sensors.sensorBiases, sensors.sensorSigmas)
print("sensors Noisy updated: ", sensors.sensorsNoisy)

sensors.update()
print("sensors True update: ", sensors.sensorsTrue)
print("sensors Noisy updated: ", sensors.sensorsNoisy)
print("updateTicks: ", sensors.updateTicks)

sensors.reset()
print("sensors True reset: ", sensors.sensorsTrue)
print("sensors Noisy reset: ", sensors.sensorsNoisy)
