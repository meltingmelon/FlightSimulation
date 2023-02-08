import math
import pickle
from ece163.Modeling import VehicleAerodynamicsModel as VAM
from ece163.Constants import VehiclePhysicalConstants as VPC
from ece163.Containers import States
from ece163.Containers import Inputs
from ece163.Containers import Controls
from ece163.Containers import Linearized
from ece163.Utilities import MatrixMath
from ece163.Utilities import Rotations

def computeGains(tuningParameters=Controls.controlTuning(), linearizedModel=Linearized.transferFunctions()):
    """
    Function to compute the control gains using the tuning parameters outline in Beard Chapter 6. Both the lateral and longitudinal gains are calculated. No check is made for frequency separation.
    Transfer function input comes from the VehiclePerturbationModels which rely on the VehicleTrim.py (provided) to compute the trim values.

    :param tuningParameters: class controlTuning from Containers.Controls
    :return controlGains: class controlGains from Containers.Controls
    """
    controlGains = Controls.controlGains()
    referenceCommands = Controls.referenceCommands()

    Vg = linearizedModel.Va_trim

    a_phi1 = linearizedModel.a_phi1
    a_phi2 = linearizedModel.a_phi2
    a_beta1 = linearizedModel.a_beta1
    a_beta2 = linearizedModel.a_beta2
    a_theta1 = linearizedModel.a_theta1
    a_theta2 = linearizedModel.a_theta2
    a_theta3 = linearizedModel.a_theta3
    a_V1 = linearizedModel.a_V1
    a_V2 = linearizedModel.a_V2
    a_V3 = linearizedModel.a_V3

    Wn_roll = tuningParameters.Wn_roll
    Zeta_roll = tuningParameters.Zeta_roll
    Wn_course = tuningParameters.Wn_course
    Zeta_course = tuningParameters.Zeta_course
    Wn_sideslip = tuningParameters.Wn_sideslip
    Zeta_sideslip = tuningParameters.Zeta_sideslip

    Wn_pitch = tuningParameters.Wn_pitch
    Zeta_pitch = tuningParameters.Zeta_pitch
    Wn_altitude = tuningParameters.Wn_altitude
    Zeta_altitude = tuningParameters.Zeta_altitude
    Wn_SpeedfromThrottle = tuningParameters.Wn_SpeedfromThrottle
    Zeta_SpeedfromThrottle = tuningParameters.Zeta_SpeedfromThrottle
    Wn_SpeedfromElevator = tuningParameters.Wn_SpeedfromElevator
    Zeta_SpeedfromElevator = tuningParameters.Zeta_SpeedfromElevator

    # Lateral Gains
    controlGains.kp_roll = ((Wn_roll ** 2)) /a_phi2
    controlGains.kd_roll = ((2 * Zeta_roll * Wn_roll) - a_phi1) / a_phi2
    controlGains.ki_roll = 0.001

    controlGains.kp_sideslip = ((2 * Zeta_sideslip * Wn_sideslip) - a_beta1) / a_beta2
    controlGains.ki_sideslip = (Wn_sideslip ** 2) / a_beta2

    controlGains.kp_course = (2 * Zeta_course * Wn_course * Vg) / VPC.g0
    controlGains.ki_course = ((Wn_course ** 2) * Vg) / VPC.g0

    #Longitudinal Gains
    # pitchErrorRatio = ((Wn_pitch ** 2) - a_theta2) / abs(a_theta3) 
    controlGains.kp_pitch = ((Wn_pitch ** 2) - a_theta2) / a_theta3 
    controlGains.kd_pitch = ((2 * Zeta_pitch * Wn_pitch) - a_theta1) / a_theta3

    kThetaDC = (controlGains.kp_pitch * a_theta3) / (a_theta2 + (controlGains.kp_pitch * a_theta3)) 
    controlGains.kp_altitude = (2 * Zeta_altitude * Wn_altitude) / (kThetaDC * Vg)
    controlGains.ki_altitude = (Wn_altitude ** 2) / (kThetaDC * Vg)

    controlGains.kp_SpeedfromThrottle = ((2 * Zeta_SpeedfromThrottle * Wn_SpeedfromThrottle) - a_V1) / a_V2
    controlGains.ki_SpeedfromThrottle = (Wn_SpeedfromThrottle ** 2) / a_V2

    controlGains.kp_SpeedfromElevator = (a_V1 - (2 * Zeta_SpeedfromElevator * Wn_SpeedfromElevator))/ (kThetaDC * VPC.g0)
    controlGains.ki_SpeedfromElevator = -1 * ((Wn_SpeedfromElevator ** 2) / (kThetaDC * VPC.g0))

    return controlGains

def computeTuningParameters(controlGains=Controls.controlGains(), linearizedModel=Linearized.transferFunctions()):
    """
    Function to compute the tuning parameters given the gains in the successive loop closure,
    needs a try block to deal with taking the square root of negative number. Function should never fail, if an exception offcures,
    return an empty(inited) tuningParameters class. Transfer function input comes from the VehiclePerturbationModels 
    which rely on the VehicleTrim.py(provided) to compute the trim values.

    :param controlGains: class controlGrains from Containers.Controls.
    :param linearizedModel: class transferFunction from Containers.Linearized

    :return tuningParameters: class controlTuning from Containers.Controls
    """

    tuningParameters = Controls.controlTuning()

    Vg = linearizedModel.Va_trim

    a_phi1 = linearizedModel.a_phi1
    a_phi2 = linearizedModel.a_phi2
    a_beta1 = linearizedModel.a_beta1
    a_beta2 = linearizedModel.a_beta2
    a_theta1 = linearizedModel.a_theta1
    a_theta2 = linearizedModel.a_theta2
    a_theta3 = linearizedModel.a_theta3
    a_V1 = linearizedModel.a_V1
    a_V2 = linearizedModel.a_V2
    a_V3 = linearizedModel.a_V3

    kp_roll = controlGains.kp_roll
    kd_roll = controlGains.kd_roll
    ki_roll = controlGains.ki_roll
    kp_sideslip = controlGains.kp_sideslip
    ki_sideslip = controlGains.ki_sideslip
    kp_course = controlGains.kp_course
    ki_course = controlGains.ki_course
    kp_pitch = controlGains.kp_pitch
    kd_pitch = controlGains.kd_pitch
    kp_altitude = controlGains.kp_altitude
    ki_altitude = controlGains.ki_altitude
    kp_SpeedfromThrottle = controlGains.kp_SpeedfromThrottle
    ki_SpeedfromThrottle = controlGains.ki_SpeedfromThrottle
    kp_SpeedfromElevator = controlGains.kp_SpeedfromElevator
    ki_SpeedfromElevator = controlGains.ki_SpeedfromElevator



    try:
        tuningParameters.Wn_roll = math.sqrt(kp_roll * a_phi2)
        Wn_roll = tuningParameters.Wn_roll
        tuningParameters.Zeta_roll = ((kd_roll * a_phi2) + a_phi1) / (2 * Wn_roll)
        
        tuningParameters.Wn_course = math.sqrt((VPC.g0 * ki_course)/ Vg )
        Wn_course = tuningParameters.Wn_course
        tuningParameters.Zeta_course = (kp_course * VPC.g0) / (2 * Wn_course * Vg)
        
        tuningParameters.Wn_sideslip = math.sqrt(a_beta2 * ki_sideslip)
        Wn_sideslip = tuningParameters.Wn_sideslip
        tuningParameters.Zeta_sideslip = (a_beta1 + (a_beta2 * kp_sideslip)) / (2 * Wn_sideslip)

        tuningParameters.Wn_pitch = math.sqrt(a_theta2 + (kp_pitch * a_theta3))
        Wn_pitch = tuningParameters.Wn_pitch
        kThetaDC = (controlGains.kp_pitch * a_theta3) / (Wn_pitch ** 2)
        tuningParameters.Zeta_pitch = (a_theta1 + (kd_pitch * a_theta3)) / (2 * Wn_pitch)
        
        tuningParameters.Wn_altitude = math.sqrt(kThetaDC * Vg * ki_altitude)
        Wn_altitude = tuningParameters.Wn_altitude
        tuningParameters.Zeta_altitude = (kThetaDC * Vg * kp_altitude) / (2 * Wn_altitude)
        
        tuningParameters.Wn_SpeedfromThrottle = math.sqrt(a_V2 * ki_SpeedfromThrottle)
        Wn_SpeedfromThrottle = tuningParameters.Wn_SpeedfromThrottle
        tuningParameters.Zeta_SpeedfromThrottle = (a_V1 + (a_V2 * kp_SpeedfromThrottle)) / (2 * Wn_SpeedfromThrottle)
        
        tuningParameters.Wn_SpeedfromElevator = math.sqrt(-1 * kThetaDC * VPC.g0 * ki_SpeedfromElevator)
        Wn_SpeedfromElevator = tuningParameters.Wn_SpeedfromElevator

        tuningParameters.Zeta_SpeedfromElevator = (a_V1 - (kThetaDC * VPC.g0 * kp_SpeedfromElevator)) / (2 * Wn_SpeedfromElevator)
        return tuningParameters

    except ValueError:
        tuningParameters = Controls.controlTuning()
        return tuningParameters

    except:
        tuningParameters = Controls.controlTuning()
        return tuningParameters


