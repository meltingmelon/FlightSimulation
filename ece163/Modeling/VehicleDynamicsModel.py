import math
from ..Containers import States
from ..Utilities import MatrixMath as mm
from ..Utilities import Rotations
from ..Constants import VehiclePhysicalConstants as VPC

# Ugh
# I left this here just in case
Gamma3 = VPC.Jzz/VPC.Jdet
Gamma4 = VPC.Jxz / VPC.Jdet
Gamma5 = (VPC.Jzz - VPC.Jxx) / VPC.Jyy
Gamma6 = VPC.Jxz / VPC.Jyy
Gamma8 = VPC.Jxx / VPC.Jdet

class VehicleDynamicsModel:
    """
    This module is where all of the vehicle dynamics are computed for the simulation. It includes the kinematics of
    both the translational and rotational dynamics. Included are both the derivative, and the integration functions, and the
    rotations of forces to the body frame.
    
    """
    def __init__(self, dT = VPC.dT):
        """
        Initialize self. See help(type(self)) for accurate signature
        """
        self. dT = dT
        self.state= States.vehicleState()
        self.dot = States.vehicleState()
        
        return

    def ForwardEuler(self, forcesMoments):
        """
        Function to do the simple forwards integration of the state using the derivative function. 
        State is integrated using the x_{k+t} =  x_{k} + dx/dt * dT.
        State is updated in place, use get VehicleState to retrive state.

        :param forceMoments: forces [N] and moments [N-m] defined in forcesMoments class
        :return: none
        """

        state = VehicleDynamicsModel.getVehicleState(self)

        #take the derivative of the state
        dot = VehicleDynamicsModel.derivative(self,state, forcesMoments)

        #integrate the state
        state = VehicleDynamicsModel.IntegrateState(self, self.dT, state, dot)

        return state

    def IntegrateState(self, dT, state, dot):
        """
        Updates the state given the derivative, and a time step. Attitude propagation is implemented as a DCM
        matrix exponential solution, all others are currently a forward integration [x]k+1 = [x]k + xdot*dT. State
        is updated in place, use getVehicleState to retrieve state.

        :param state: vehicle state
        :param dot: derivative of vehicle state
        :return: newState()
        """
        # NOTE: This implementation follows eqns from Attitude and Kinematic cheat sheets


        # get Rexp
        R = VehicleDynamicsModel.Rexp(self, dT,state,dot)

        # matrix multiply to get rotation matrix
        newR = mm.matrixMultiply(R, state.R)
        
       # set's new newState's R to the R value calculated one
        # newState.R = newR
        R = newR

        # update euler angles using DCM2euler
        # newState.yaw, newState.pitch, newState.roll = Rotations.dcm2Euler(newState.R)
        yaw, pitch, roll = Rotations.dcm2Euler(R)

        # manual forward integration using eqn for all other  variables
        # newState.pn = state.pn + dot.pn * dT
        # newState.pe = state.pe + dot.pe * dT
        # newState.pd = state.pd + dot.pd * dT

        pn = state.pn + (dot.pn * dT)
        pe = state.pe + (dot.pe * dT)
        pd = state.pd + (dot.pd * dT)

        u = state.u + (dot.u * dT)
        v = state.v + (dot.v  * dT)
        w = state.w + (dot.w * dT)

        p = state.p + (dot.p * dT)
        q = state.q + (dot.q * dT)
        r = state.r + (dot.r * dT)

        newState = States.vehicleState(pn,pe,pd,u,v,w,yaw,pitch,roll,p,q,r,R)
        newState.Va = state.Va
        newState.alpha = state.alpha
        newState.beta = state.beta
        newState.chi = state.chi
        return newState

    def Rexp(self, dT,state,dot):
        """
        Performs the matrix exponential closed form solution for the DCM integration from body-fixed rates.
        Uses the Murray formulation for the closed form solution.
        Implements trapezoidal rate for integration.

        :param state: vehicle state
        :param dot: derivative of vehicle state
        :param dT: Time step[s]
        :return Rexp: the matrix exponential to update the state
        """
        
        # Create Identity Matrix [3x3]
        identityMatrix  = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]]

        #Skew symmetric matrix [3x3 ] 
        #  [ 0  -r     q ]
        #  [  r   0   -p ]
        #  [-q  p     0 ] 

        angRateDerive = [[dot.p], [dot.q], [dot.r]]
        currAngRates = [[state.p], [state.q], [state.r]]

        angRateDisplace = mm.matrixScalarMultiply((dT / 2), angRateDerive)
        newAngRates = mm.matrixAdd(currAngRates, angRateDisplace)
        # print(newAngRates)
        
        newP = newAngRates[0][0]
        newQ = newAngRates[1][0]
        newR = newAngRates[2][0]

        skewSym = mm.matrixSkew(newP, newQ, newR )

        # Distance of skew symmetric matrix
        skewMag = math.hypot(newP, newQ, newR)

        # Skew coefficients
        if skewMag <= 0.2:
            skewCoeff = dT - (((dT ** 3) * skewMag ** 2) / 6) + ((dT ** 5 * skewMag ** 4) / 120)
            skewOffCoeff = ((dT ** 2) / 2) - ((dT ** 4 * skewMag ** 2) / 24) + ((dT ** 6 * skewMag ** 4) / 720)
        else:
            skewCoeff = (math.sin(skewMag * dT)) / skewMag
            skewOffCoeff = (1 - math.cos(skewMag * dT)) / (skewMag ** 2)
        
        # Split into components
        skewSquared = mm.matrixMultiply(skewSym,skewSym)
        skew1 = mm.matrixScalarMultiply(skewCoeff, skewSym)
        skew2 = mm.matrixScalarMultiply(skewOffCoeff, skewSquared)
        skewDiff = mm.matrixSubtract(identityMatrix, skew1) 

        # Find Rexp
        Rexp = mm.matrixAdd(skewDiff, skew2)

        return Rexp

    def Update(self,forcesMoments):
        """
        Function that implements the integration such that the state is updated using the forces and moments passed in as the arguments, and the time step dT.
        State is updated in place (self.state is updated). Use getVehicleState to retrive state.
        Time step is defined in VehiclePhysicalConstants.py

        :param forcesMoments: forces [N] and moments [N-m] defined in forcesMoments class
        :return: none
        """
        # Simply updates the internal state using ForwardEuler
        self.state = VehicleDynamicsModel.ForwardEuler(self, forcesMoments)

        return
    
    def derivative(self,state,forcesMoments):
        """ 
            Function to compute the derivative of the state given body frame forces and moments.

            :param state: the vehicle state
            :param forcesMoments: forces [N] and moments [N-m] defined in forcesMoments class
            :returns dot: the state derivative
        """

        # Derivative  eqns found in Chapter 3 Summary
        # and Attitude Cheat Sheet for Derivative of 

        # derivative of positions
        inertialVelocity = [[state.u],[state.v],[state.w]] 
        transposeR = mm.matrixTranspose(state.R)

        translationVelocity  = mm.matrixMultiply(transposeR, inertialVelocity)

        # dot.pn = translationVelocity[0][0]
        # dot.pe = translationVelocity[1][0]
        # dot.pd = translationVelocity[2][0]

        pn = translationVelocity[0][0]
        pe = translationVelocity[1][0]
        pd = translationVelocity[2][0]
        # derivative of euler angles
        angRates = [[state.p],[state.q], [state.r]]

        R = [
            [1, math.sin(state.roll) * math.tan(state.pitch), math.cos(state.roll) * math.tan(state.pitch)],
            [0, math.cos(state.roll), -1 * math.sin(state.roll)],
            [0, math.sin(state.roll) / math.cos(state.pitch), math.cos(state.roll) / math.cos(state.pitch)]]
        
        eulerDerivatives = mm.matrixMultiply(R,angRates)
        skewSym = mm.matrixSkew(state.p, state.q, state.r)

        # dot.roll = eulerDerivatives[0][0]
        # dot.pitch = eulerDerivatives[1][0]
        # dot.yaw = eulerDerivatives[2][0]

        roll = eulerDerivatives[0][0]
        pitch = eulerDerivatives[1][0]
        yaw = eulerDerivatives[2][0]

        # derivative R
        # dot.R = mm.matrixScalarMultiply(-1, (mm.matrixMultiply(skewSym, state.R)))
        R = mm.matrixScalarMultiply(-1, (mm.matrixMultiply(skewSym, state.R)))
        
        # derivative of velocity vector
        fxM = forcesMoments.Fx / VPC.mass
        fyM = forcesMoments.Fy / VPC.mass
        fzM = forcesMoments.Fz / VPC.mass

        # velocity vector
        f = [
            [fxM],
            [fyM],
            [fzM]]
        
        displacements = [
            [(state.r * state.v) - (state.q * state.w)],
            [(state.p * state.w) - (state.r * state.u)],
            [(state.q * state.u) - (state.p * state.v)]]

        velocityDerivatives = mm.matrixAdd(displacements, f)

        # dot.u = velocityDerivatives[0][0]
        # dot.v = velocityDerivatives[1][0]
        # dot.w = velocityDerivatives[2][0]

        u = velocityDerivatives[0][0]
        v = velocityDerivatives[1][0]
        w = velocityDerivatives[2][0]
        
        # derivative of angular rates
        angDisplacement = [
            [(VPC.Gamma1 * state.p * state.q) - (VPC.Gamma2 * state.q * state.r)],
            [(Gamma5 * state.p * state.r) - (Gamma6 * (state.p ** 2 - state.r ** 2))],
            [(VPC.Gamma7 * state.p * state.q) - (VPC.Gamma1 * state.q * state.r)]]
        
        angOffset = [
            [(Gamma3 * forcesMoments.Mx) + (Gamma4 * forcesMoments.Mz)],
            [forcesMoments.My / VPC.Jyy],
            [(Gamma4 * forcesMoments.Mx) + (Gamma8 * forcesMoments.Mz) ]]
        
        angRateDerivative = mm.matrixAdd(angDisplacement, angOffset)

        # dot.p = angRateDerivative[0][0]
        # dot.q = angRateDerivative[1][0]
        # dot.r = angRateDerivative[2][0]

        p = angRateDerivative[0][0]
        q = angRateDerivative[1][0]
        r = angRateDerivative[2][0]        

        dot = States.vehicleState(pn,pe,pd,u,v,w,yaw,pitch,roll,p,q,r)
        dot.R = R
        return dot

    def getVehicleState(self):
        """
        Wrapper function to return the vehicle state
        :returns  state:  from class vehicleState
        """
        return self.state

    def reset(self):
        """
        Resets module to original states so it can run again
        :returns: internal states are reset
        """
        self.dot = States.vehicleState()
        self.state = States.vehicleState()
        return

    def resetVehicleState(self):
        """
        Wrapper function to reset the vehicle state to initial conditions
        :returns: internal state is reset to initial conditions
        """
        self.state= States.vehicleState()
        return

    def setVehicleState(self,state):
        """
        Wrapper function to set the vehicle state
        :param state: state to be set
        :returns: internal state is set to match arguments
        """
        self.state= state
        return
