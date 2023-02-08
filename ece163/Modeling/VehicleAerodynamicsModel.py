import math
from ..Containers import States
from ..Containers import Inputs
from ..Modeling import VehicleDynamicsModel as VD
from ..Modeling import WindModel
from ..Utilities import MatrixMath as mm
from ..Constants import VehiclePhysicalConstants as VPC

class VehicleAerodynamicsModel:
    """
    Vehicle Aerodynamics Model will take the parameters from the VehiclePhysicalConstants and use the vehicle state to
    compute the forces and moments about the vehicle for integration. These are linearized aerodynamics using stability
    derivatives and simplified aerodynamics; the only deviation from this is the post-stall lift model used to convert lift
    from the linear model to a flat plate at an angle.

    Note that stability derivatives are non-dimensional, and thus need to be converted back to actual forces and moments

    """
    def __init__(self, initialSpeed=VPC.InitialSpeed, initialHeight=VPC.InitialDownPosition):
        """
        Initialization of the internal classes which are used to track the vehicle aerodynamics and dynamics.
        """
        self.vehicleDynamics = VD.VehicleDynamicsModel()
        self.vehicleDynamics.state.u = initialSpeed
        self.vehicleDynamics.state.pd = initialHeight
        self.windModel = WindModel.WindModel()
        # self.wind = States.windState()
        return

    def CalculateAirspeed(self, state, wind):
        """
        Calculates the total airspeed, as well as angle of attack and side-slip angles from the wind and current state.
        Needed for further aerodynamic force calculations. Va, wind speed [m/s], alpha, angle of attack [rad], and
        beta, side-slip angle [rad] are updated in the state class.

        :param state: current vehicle state (need the velocities)
        :param wind: current wind state (global and gust)
        :return: none
        """

        #Uses equations in CH4.4 in the book

        Wg = [[wind.Wu], [wind.Wv], [wind.Ww]]
        Ww =[[wind.Wn], [wind.We], [wind.Wd]]
        Vg = [[state.u], [state.v], [state.w]]

        Ws = math.hypot(wind.Wn, wind.We, wind.Wd)
        if (Ws == 0.0):
            Gammaw  = 0.0
        else:
            Gammaw = -math.asin(wind.Wd / Ws )
        
        Xw = math.atan2(wind.We, wind.Wn)

        RAzimuthElevation = [
            [math.cos(Xw) * math.cos(Gammaw), math.sin(Xw) * math.cos(Gammaw), -1 * math.sin(Gammaw)], 
            [-1 * math.sin(Xw), math.cos(Xw), 0],
            [math.cos(Xw) * math.sin(Gammaw), math.sin(Xw) * math.sin(Gammaw), math.cos(Gammaw)]]
        
        iWg = mm.matrixMultiply(mm.matrixTranspose(RAzimuthElevation), Wg)
        Wsum = mm.matrixAdd(Ww, iWg)
        VwmatrixMult = mm.matrixMultiply(state.R, Wsum)

        VaB = mm.matrixSubtract(Vg, VwmatrixMult)

        uR = VaB[0][0]
        vR = VaB[1][0]
        wR = VaB[2][0]

        Va = math.hypot(uR, vR, wR)
        alpha = math.atan2(wR, uR)

        if math.isclose(Va , 0.0):
            beta = 0.0
        else:
            beta = math.asin(vR/Va)
        
        return Va, alpha, beta

    def CalculateCoeff_alpha(self, alpha):
        """
        Function to calculate the Coefficient of Lift and Drag as a function of angle of attack. Angle of attack
        (alpha) in [rad] is contained within the state.alpha and calculated within the CalculateAirspeed function.

        :return: CL_alpha, CD_alpha, CM_alpha (all unitless)
        """


        # Calculating Sigmoid Blending Function. (eqn 4.10)
        num =( 1 + math.exp(-VPC.M * (alpha - VPC.alpha0))) + (math.exp(VPC.M * (alpha + VPC.alpha0)))
        den = (1 + (math.exp(-VPC.M * (alpha - VPC.alpha0)))) * (1 + (math.exp(VPC.M * (alpha + VPC.alpha0))))

        sigmoid= num / den

        # Calculating Coefficients that account for beyond stall. Referenced from Gabe's lecture notes.
        CL= VPC.CL0 + (VPC.CLalpha * alpha)
        CL_alpha = ((1 - sigmoid)  * CL) + (sigmoid * (2 * math.sin(alpha) * math.cos(alpha)))

        CD = VPC.CDp + (((CL_alpha * alpha) ** 2) / (math.pi  * VPC.AR * VPC.e))
        CD_alpha = ((1 - sigmoid) * CD) + (sigmoid * (2 * ((math.sin(alpha) ** 2))))
        
        CM_alpha = VPC.CM0 + (VPC.CMalpha * alpha)

        return CL_alpha, CD_alpha, CM_alpha

    def CalculatePropForces(self,Va,throttle):
        """
        Function to calculate the propeller forces and torques on the
        aircraft. Requires the airspeed (Va) in [m/s] in state.Va to perform the calculations. Uses the fancy propeller
        model that parameterizes the torque and thrust coefficients of the propeller using the advance ratio. See
        ECE163_PropellerCheatSheet.pdf for details.

        :param Throttle: Throttle Input [0-1]
        :return Fx_prop [N], Mx_prop [N-m] 
        """  

        
        # This was all derived from the equations in the Propeller Cheat Sheet

        KT = KE  = 60 / (2 * math.pi * VPC.KV)
        vin = VPC.V_max * throttle

        # I split the variables for omega to make my  life less hellish
        a = (VPC.rho * (VPC.D_prop ** 5) * VPC.C_Q0) / (4 * (math.pi ** 2))
        b = ((VPC.rho * (VPC.D_prop ** 4) * Va  * VPC.C_Q1) / (2 * math.pi)) + ((KT * KE) / VPC.R_motor)
        c = (VPC.rho * (VPC.D_prop ** 3) * (Va ** 2) * VPC.C_Q2) - (KT * ((vin) / VPC.R_motor) ) +  (KT * VPC.i0)
        
        if  ((b ** 2)  < (4 * a * c)):
            omega = 100
        else:
            omega = ((-1 * b) + math.sqrt((b ** 2) - (4 * a * c))) / (2 * a)


        J = (2 * math.pi * Va) / (omega * VPC.D_prop) 

        cT = VPC.C_T0 + (VPC.C_T1 * J) + (VPC.C_T2 * (J ** 2))
        cQ = VPC.C_Q0 + (VPC.C_Q1 * J) + (VPC.C_Q2 * (J ** 2))

        Fx_prop = (VPC.rho * (omega ** 2) * (VPC.D_prop ** 4) * cT) / (4 * (math.pi ** 2))
        Mx_prop = (-1)* ((VPC.rho * (omega ** 2) * (VPC.D_prop ** 5) * cQ) / (4 * (math.pi ** 2)))

        return Fx_prop, Mx_prop

    def Update(self, controls):
        """
        Function that uses the current state (internal), wind (internal), and controls (inputs) to calculate the forces,
        and then do the integration of the full 6-DOF non-linear equations of motion. Wraps the VehicleDynamic-
        sModel class as well as the windState internally. The Wind and the vehicleState are maintained internally.

        :param controls: controlInputs class (Throttle, Elevator, Aileron, Rudder)
        :return none: state is updated internally
        """

        state = VehicleAerodynamicsModel.getVehicleState(self)
        wind = VehicleAerodynamicsModel.getWindState(self)

        newForces = VehicleAerodynamicsModel.updateForces(self, state, wind, controls)
        self.vehicleDynamics.Update(newForces)
        return

    def aeroForces(self, state):
        """
        Function to calculate the Aerodynamic Forces and Moments using the lin-
        earized simplified force model and the stability derivatives in VehiclePhysicalConstants.py file. Specifi-
        cally does not include forces due to control surface deflection. Requires airspeed (Va) in [m/s], angle of
        attack (alpha) in [rad] and sideslip angle (beta) in [rad] from the state.

        :param state: current vehicle state (need the velocites)
        :return: Aerodynamic forces (forcesMoments class)
        """
        # aeroForces = Inputs.forcesMoments()

        # A conditional when there is no airspeed
        # if (state.Va == 0):
        #     return aeroForces
        if (state.Va == 0):
            aeroForces.Fx = 0.0
            aeroForces.Fy = 0.0
            aeroForces.Fz = 0.0
            aeroForces.Mx = 0.0
            aeroForces.My = 0.0
            aeroForces.Mz = 0.0

        else:
            alphaCL, alphaCD, alphaCM = VehicleAerodynamicsModel.CalculateCoeff_alpha(self, state.alpha) 

            R = [
                [math.cos(state.alpha), -1 * math.sin(state.alpha)],
                [math.sin(state.alpha), math.cos(state.alpha)]
                ]

            qBar = (VPC.c * state.q) / (2 * state.Va)
            pBar = (VPC.b  * state.p) / (2 * state.Va)
            rBar = (VPC.b * state.r) / (2 * state.Va)

            pitchRateCL = VPC.CLq * qBar
            pitchRateCD = VPC.CDq * qBar
            
            k = (1 / 2) * VPC.rho * (state.Va ** 2) * VPC.S

            fLift = k * (alphaCL+ pitchRateCL)
            fDrag = k * (alphaCD + pitchRateCD)

            fLD = [[-1 * fDrag], [-1 * fLift]]
            fLateral = mm.matrixMultiply(R,fLD)

            # Aero Forces
            Fx = fLateral[0][0]
            Fy = k * (VPC.CY0 + (VPC.CYbeta * state.beta) + (VPC.CYp * pBar) + (VPC.CYr * rBar) )
            Fz = fLateral[1][0]

            # Aero Moments
            Mx = (k * VPC.b) * (VPC.Cl0 + (VPC.Clbeta * state.beta) + (VPC.Clp * pBar) + (VPC.Clr * rBar))
            My = (k * VPC.c) * ((VPC.CM0 + (VPC.CMalpha * state.alpha))+ (VPC.CMq * qBar))
            Mz = (k * VPC.b) * (VPC.Cn0 + (VPC.Cnbeta * state.beta) + (VPC.Cnp * pBar) + (VPC.Cnr * rBar))


            aeroForces = Inputs.forcesMoments(Fx,Fy,Fz,Mx,My,Mz)

        return aeroForces

    def controlForces(self, state, controls):
        """
        Function to calculate aerodynamic forces from control surface deflec-
        tions (including throttle) using the linearized aerodynamics and simplified thrust model. Requires airspeed
        (Va) in [m/s] and angle of attack (alpha) in [rad] both from state.Va and state.alpha respectively.

        :param controls: inputs to the aircraft controlInputs()
        :return: control surface forces forcesMoments class
        """
        
        controlForce =Inputs.forcesMoments()
        alpha = state.alpha

        R = [
            [math.cos(state.alpha), -math.sin(state.alpha) ],
             [math.sin(state.alpha), math.cos(state.alpha)]]

        k = (1 / 2) * VPC.rho * (state.Va ** 2) * VPC.S

        fLiftC = k * (VPC.CLdeltaE * controls.Elevator)
        fDragC = k *  (VPC.CDdeltaE * controls.Elevator)

        fLD = [[-1 * fDragC], [-1 * fLiftC]]
        fLateralC = mm.matrixMultiply(R,fLD)

        # control Forces

        controlFx = fLateralC[0][0]
        controlFy = k * ((VPC.CYdeltaA * controls.Aileron) + (VPC.CYdeltaR * controls.Rudder)) 
        controlFz = fLateralC[1][0]


        # control Moments
        controlMx = (k * VPC.b) * ((VPC.CldeltaA * controls.Aileron) + (VPC.CldeltaR * controls.Rudder))
        controlMy = (k * VPC.c) * (VPC.CMdeltaE * controls.Elevator)
        controlMz = (k * VPC.b) * ((VPC.CndeltaA * controls.Aileron) + (VPC.CndeltaR * controls.Rudder))

        propForces, propMoments = VehicleAerodynamicsModel.CalculatePropForces(self, state.Va, controls.Throttle)

        # controlForce.Fx = controlFx + propForces
        # controlForce.Fy = controlFy 
        # controlForce.Fz = controlFz

        # controlForce.Mx = controlMx + propMoments
        # controlForce.My = controlMy
        # controlForce.Mz = controlMz

        Fx = controlFx + propForces
        Fy = controlFy 
        Fz = controlFz

        Mx = controlMx + propMoments
        My = controlMy
        Mz = controlMz

        controlForce =Inputs.forcesMoments(Fx, Fy, Fz, Mx, My, Mz)

        return controlForce

    def getVehicleDerivative(self):
        """
        Wrapper function to return vehicle state derivative from module

        :return derivative: class of vehicleState
        """

        return self.vehicleDynamics.dot

    def getVehicleState(self):
        """
        Wrapper function to return vehicle state form module
        :return vehicle state class:
        """
        return self.vehicleDynamics.state

    def getWindState(self):
        """
        Wrapper function to return wind state from module
        :return wind state class:
        """
        return self.windModel.Wind

    def setVehicleDerivative(self, dot):
        """
        Wrapper function to set the vehicle state derivative from outside module

        :param dot: class of vehicleState
        :return none:
        """
        self.vehicleDynamics.dot = dot
        return

    def gravityForces(self, state):
        """
        Function to project gravity forces into the body frame. Uses the gravity
        constant g0 from physical constants and the vehicle mass. Fg = m * R * [0 0 g0]’

        :param state: current vehicle state (need  the rotation matrix)
        :return gravity forces: forcesMoments class
        """

        gravityVector = [[0], [0], [VPC.mass * VPC.g0]]

        # Try multiplying scalar mass with rotation matrix first

        fG = mm.matrixMultiply(state.R, gravityVector)
        # gravity.Fx = fG[0][0]
        # gravity.Fy = fG[1][0]
        # gravity.Fz = fG[2][0]

        Fx = fG[0][0]
        Fy = fG[1][0]
        Fz = fG[2][0]
        gravity = Inputs.forcesMoments(Fx, Fy, Fz)

        return gravity

    def reset(self):
        """
        Resets module to its original state so it can run again

        :return none:
        """
        self.vehicleDynamics = VD.VehicleDynamicsModel()
        self.vehicleDynamics.state.u = VPC.InitialSpeed
        self.vehicleDynamics.state.pd = VPC.InitialDownPosition
        self.windModel = WindModel.WindModel()
        return

    def setVehicleState(self, state):
        """
        Wrapper function to set the vehicle state from outside module

        :param state: class of vehicleState
        :return: none
        """
        self.vehicleDynamics.state = state
        return

    def setWindModel(self, Wn=0.0, We=0.0, Wd=0.0, drydenParamters= Inputs.drydenParameters(Lu=0.0, Lv=0.0, Lw=0.0, sigmau=0.0, sigmav=0.0, sigmaw=0.0)):
        """
        Wrapper function that will inject constant winds and gust parameters into the wind model using the con-
        stant wind in the inertial frame (steady wind) and gusts that are stochastically derived in the body frame
        using the Dryden wind gust models.

        :return: none
        """
        self.windModel.Wind.Wn = Wn
        self.windModel.Wind.We = We
        self.windModel.Wind.Wd = Wd
        self.windModel.drydenParamters = drydenParamters
        
        return

    def updateForces(self, state, wind, controls):
        """
        Function to update all of the aerodynamic, propulsive, and
        gravity forces and moment. Self is updated with new values using the class definition

        :param state: current vehicle state
        :param wind: current environmental wind
        :param controls: current vehicle control surface deflections
        :return: total forces (forcesMoments class)
        """

        # totalForces = Inputs.forcesMoments()

        Va, alpha, beta = self.CalculateAirspeed(state,wind)
        state.Va = Va
        state.alpha = alpha
        state.beta = beta

        gravity = VehicleAerodynamicsModel.gravityForces(self, state)
        aero = VehicleAerodynamicsModel.aeroForces(self, state)
        control = VehicleAerodynamicsModel.controlForces(self, state, controls)

        # totalForces.Fx = gravity.Fx +aero.Fx + control.Fx
        # totalForces.Fy = gravity.Fy + aero.Fy + control.Fy
        # totalForces.Fz = gravity.Fz + aero.Fz + control.Fz
        # totalForces.Mx = gravity.Mx + aero.Mx + control.Mx
        # totalForces.My = gravity.My + aero.My + control.My
        # totalForces.Mz = gravity.Mz + aero.Mz + control.Mz

        Fx = gravity.Fx +aero.Fx + control.Fx
        Fy = gravity.Fy + aero.Fy + control.Fy
        Fz = gravity.Fz + aero.Fz + control.Fz
        Mx = gravity.Mx + aero.Mx + control.Mx
        My = gravity.My + aero.My + control.My
        Mz = gravity.Mz + aero.Mz + control.Mz

        totalForces = Inputs.forcesMoments(Fx, Fy, Fz, Mx, My, Mz)

        return totalForces

