import math
import pickle
from ece163.Modeling import VehicleAerodynamicsModel
from ece163.Constants import VehiclePhysicalConstants as VPC
from ece163.Containers import States
from ece163.Containers import Inputs
from ece163.Containers import Linearized
from ece163.Utilities import MatrixMath
from ece163.Utilities import Rotations
from ece163.Controls import VehicleTrim


def CreateTransferFunction(trimState, trimInputs):
    """
    Function to fill the transfer function parameters used for the successive loop closure from the given trim
    state and trim inputs. Note that these parameters will be later used to generate actual control loops. Pa-
    rameters are updated in self.transferFunction

    :param trimState: vehicle trim state, as calculated in VehicleTrim code
    :param trimInputs:  vehicle trim inputs, as calculated in VehicleTrim code

    :return: none, parameters are updated internally in self.transferFunction 
    """
    
    transferFunction = Linearized.transferFunctions()

    # set mass
    m = VPC.mass
    Va = trimState.Va
    alphaTrim = trimState.alpha
    pitchTrim = trimState.pitch
    throttleTrim = trimInputs.Throttle
    elevatorTrim = trimInputs.Elevator

    a_phi1 =  ((-1 / 2) * VPC.rho * (Va ** 2) * VPC.S * VPC.b * VPC.Cpp) * (VPC.b / (2 * Va))
    a_phi2 = (1 / 2) * VPC.rho * (Va ** 2) * VPC.S * VPC.b * VPC.CpdeltaA
    
    a_beta1 = -1 * (((VPC.rho * Va * VPC.S) / (2 * m)) * VPC.CYbeta)
    a_beta2 = (((VPC.rho * Va * VPC.S) / (2 * m)) * VPC.CYdeltaR)

    a_theta1 = (-1)*((VPC.rho*(Va**2)*VPC.S*VPC.c)/(2*VPC.Jyy))*VPC.CMq*(VPC.c/(2*Va)) 
    a_theta2 = -1 * (((VPC.rho * (Va ** 2) * VPC.c * VPC.S) / (2 * VPC.Jyy)) * VPC.CMalpha)
    a_theta3 = ((VPC.rho * (Va ** 2) * VPC.c * VPC.S) / (2 * VPC.Jyy)) * VPC.CMdeltaE

    dTdV = dThrust_dVa(Va,throttleTrim, 0.5)
    dTdT = dThrust_dThrottle(Va, throttleTrim, 0.01)
    
    a_V1 = (((VPC.rho*Va*VPC.S)/m)*(VPC.CD0+(VPC.CDalpha*alphaTrim)+(VPC.CDdeltaE*elevatorTrim)))-((1/m)*dTdV) 
    a_V2 = (1 / m) * dTdT
    a_V3 = VPC.g0 * (math.cos(pitchTrim - alphaTrim))

    transferFunction.a_phi1 = a_phi1
    transferFunction.a_phi2 = a_phi2

    transferFunction.a_beta1 = a_beta1
    transferFunction.a_beta2 = a_beta2

    transferFunction.a_theta1 = a_theta1
    transferFunction.a_theta2 = a_theta2
    transferFunction.a_theta3 = a_theta3

    transferFunction.a_V1 = a_V1
    transferFunction.a_V2 = a_V2
    transferFunction.a_V3 = a_V3

    transferFunction.Va_trim = trimState.Va
    transferFunction.alpha_trim = trimState.alpha
    
    if (trimState.Va == 0):
        transferFunction.beta_trim = math.copysign((math.pi / 2), trimState.u)
    else:   
        transferFunction.beta_trim = trimState.beta

    transferFunction.gamma_trim = (trimState.pitch - trimState.alpha)
    transferFunction.theta_trim = trimState.pitch
    transferFunction.phi_trim = trimState.roll

    return  transferFunction

def dThrust_dThrottle(Va, Throttle, epsilon=0.01):
    """
    Function to calculate the numerical partial derivative of pro-
    peller thrust to change in throttle setting using the actual prop function from complex propeller model.
    
    :param epsilon: step to take for perturbation in throttle setting
    :return: dTdDelta: partial derivative [N/PWM]
    """
    VAD = VehicleAerodynamicsModel.VehicleAerodynamicsModel()

    FxProp ,  MxProp = VAD.CalculatePropForces(Va, Throttle)
    FxPropE, MxPropE = VAD.CalculatePropForces(Va, (Throttle + epsilon))

    dTdDeltaT = (FxPropE - FxProp) / epsilon

    return dTdDeltaT

def dThrust_dVa(Va, Throttle, epsilon = 0.5):
    """
    Function to calculate the numerical partial derivative of propeller
    thrust to change in airspeed using the actual prop function from complex propeller model.
    
    :param epsilon: step to take for perturbation in Velocity
    :return: dTdVa: partial derivative [N-s/m]
    """
    VAD = VehicleAerodynamicsModel.VehicleAerodynamicsModel()

    FxProp, MxProp = VAD.CalculatePropForces(Va, Throttle)
    FxPropE, MxPropE = VAD.CalculatePropForces((Va + epsilon), Throttle)

    dTdVa = (FxPropE - FxProp) / epsilon

    return dTdVa
    
    