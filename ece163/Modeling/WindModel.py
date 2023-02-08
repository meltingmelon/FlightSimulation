import math
import random
from ..Containers import Inputs
from ..Containers import States
from ..Utilities import MatrixMath as mm
from ..Constants import VehiclePhysicalConstants as VPC

class WindModel:
    """
    WindModel.py implements the Dryden gust model and allows for an update of the wind at every time step when driven
    by white noise. The constants for the wind model are embedded in the init function, and an init flag tells the __init__
    function which set of paramters to grab.
    All of the constants that are used in the Dryden Wind Gust Models. These are detailed on Beard Chapter 4.4. The
    Dryden models are typically implemented assuming a constant nominal airspeed (taken as initial speed from physical
    constants). The parameters for the Dryden gust model are defined in MIL-F-8785C.
    """
    def __init__(self, dT = 0.01, Va = 25.0, drydenParamters = VPC.DrydenNoWind):
        """
        function to initialize the wind model code. Will load the appropriate constants that parameterize the wind
        gusts from the Dryden gust model. Creates the discrete transfer functions for the gust models that are used
        to update the local wind gusts in the body frame. These are added to the inertial wind (Wn, We, Wd) that
        are simply constants. Discrete models are held in self and used in the Update function.
        To turn off gusts completely, use the DrydenNoGusts parameters.

        :param dT: time step [s] for numerical integration of wind
        :param Va: nominal flight speed [m/s]
        :drydenParamters: class from Inputs
        """
        
        self.drydenParamters = drydenParamters
        self.dT = dT
        self.Va = Va
        self.Wind = States.windState()

        self.Phi_u= [[0.0]]
        self.Gamma_u= [[0.0]]
        self.H_u = [[0.0]]

        self.Phi_v = [[0.0, 0.0], [0.0, 0.0]]
        self.Gamma_v= [[0.0],[0.0]]
        self.H_v = [[0.0,0.0]]

        self.Phi_w = [[0,0], [0, 0]]
        self.Gamma_w= [[0], [0]]
        self.H_w =[[0,0]]

        WindModel.CreateDrydenTransferFns(self, dT, Va, drydenParamters)

        self.x_u= [[0]]
        self.x_v = [[0.0], [0.0]]
        self.x_w = [[0.0], [0.0]]

        return

    def CreateDrydenTransferFns(self, dT, Va, drydenParamters):
        """
        Function creates the Dryden transfer functions in discrete form. These are used in generating the gust
        models for wind gusts (in body frame).

        :param dT: time step [s]
        :param Va: nominal flight speed [m/s]
        :param drydenParamters: Dryden Wing Gust Model from VehiclePhysicalConstants

        :return: none
        """

        if (drydenParamters == VPC.DrydenNoWind):
            self.Phi_u =[[1.0]]
            self.H_u = [[1.0]]
            self.Phi_v = [[1.0, 0.0], [0.0, 1.0]]
            self.H_v = [[1.0, 1.0]]
            self.Phi_w = [[1.0, 0.0], [0.0, 1.0]]
            self.H_w = [[1.0, 1.0]]

        
        else:
            Lu = drydenParamters.Lu
            Lv = drydenParamters.Lv
            Lw = drydenParamters.Lw

            sigmaU =drydenParamters.sigmau
            sigmaV = drydenParamters.sigmav
            sigmaW = drydenParamters.sigmaw 

            # Hu(s)
            self.Phi_u = [[math.exp( -1 * (Va / Lu) * dT)]]
            self.Gamma_u = [[(Lu / Va) * (1 - (math.exp(-1 * (Va / Lu) * dT)))]]
            self.H_u = [[sigmaU * math.sqrt((2 * Va) / (math.pi * Lu))]]

            #Hv(s)
            # Matrices
            phiMatrixV = [
                [1- ((Va / Lv) * dT), -1 * ((Va/Lv) ** 2) * dT],
                [dT, 1 + ((Va / Lv) * dT)]]
            gammaMatrixV = [[dT], [(((Lv / Va) ** 2) * (math.exp((Va / Lv) *dT) - 1)) - ((Lv / Va) * dT)]]
            hMatrixV = [[1, Va / (math.sqrt(3) * Lv)]]
            hSqrtV = sigmaV * math.sqrt((3 * Va) / (math.pi * Lv))

            expV = math.exp(-1 * (Va / Lv)* dT)

            #Hv(s) discrete values
            self.Phi_v = mm.matrixScalarMultiply(expV, phiMatrixV)
            self.Gamma_v = mm.matrixScalarMultiply(expV, gammaMatrixV)
            self.H_v = mm.matrixScalarMultiply(hSqrtV, hMatrixV)

            #Hw(s)
            # Matrices
            phiMatrixW= [
                [1- ((Va / Lw) * dT), -1 * ((Va/Lw) ** 2) * dT],
                [dT, 1 + ((Va / Lw) * dT)]]
            gammaMatrixW = [[dT], [(((Lw / Va) ** 2) * (math.exp((Va / Lw) *dT) - 1)) - ((Lw / Va) * dT)]]
            hMatrixW = [[1, Va / (math.sqrt(3) *  Lw)]]
            expW = math.exp(-1 * (Va / Lw)* dT)

            # Hw(s) discrete values
            self.Phi_w = mm.matrixScalarMultiply(expW, phiMatrixW)
            self.Gamma_w = mm.matrixScalarMultiply(expW, gammaMatrixW)
            self.H_w = mm.matrixScalarMultiply((sigmaW * (math.sqrt((3 * Va) / (math.pi * Lw)))), hMatrixW)

        return

    def getDrydenTransferFns(self):
            """
            Wrapper function to return the internals of the Dryden Transfer function in order to be able to test the code without requiring consistent internal names. 
            Returns the discretized version of the Drydem gust model as outlined in the ECE163_DrydenWindModel handout (Phi_u, Gamma_u, H_u, Phi_v, Gamma_v, H_v, Phi_w, Gamma_w, H_w)

            :params: none
            :return:  Phi_u, Gamma_u, H_u, Phi_v, Gamma_v, H_v, Phi_w, Gamma_w, H_w
            """

            return self.Phi_u, self.Gamma_u, self.H_u, self.Phi_v, self.Gamma_v, self.H_v, self.Phi_w, self.Gamma_w, self.H_w

    def Update(self, uu= None, uv=None, uw=None):
        """
        Function that updates the wind gusts and inserts them back into the .Wind portion of self. This is done by
        running white noise [Gaussian(0,1)] through the coloring filters of the Dryden Wind Gust model.
        
        :param  uu: optional argument for injecting input to Hu(s)
        :param uv: optional argument for injecting input to Hv(s), defaults to None
        :param uw: – optional argument for injecting input to Hw(s)
        
        :return: none
        """

        mu  = 0.0
        sigma = 1.0

        if  (uu == None):
            uu = random.gauss(mu,sigma)
        if (uv == None):
            uv = random.gauss(mu, sigma)
        if (uw == None):
            uw = random.random(mu, sigma)
        
        xU = self.x_u
        xV = self.x_v
        xW = self.x_w
        
        # updates Wu
        newXU = mm.matrixAdd(mm.matrixMultiply(self.Phi_u, xU), mm.matrixScalarMultiply(uu, self.Gamma_u))
        Wu = mm.matrixMultiply(self.H_u, newXU)
        self.Wind.Wu = Wu[0][0]
        self.x_u= newXU

        # updates Wv
        newXV = mm.matrixAdd(mm.matrixMultiply(self.Phi_v, xV),mm.matrixScalarMultiply(uv, self.Gamma_v))
        Wv = mm.matrixMultiply(self.H_v, newXV)
        self.Wind.Wv = Wv[0][0]
        self.x_v = newXV

        # updates Ww
        newXW = mm.matrixAdd(mm.matrixMultiply(self.Phi_w, xW), mm.matrixScalarMultiply(uw, self.Gamma_w))
        Ww = mm.matrixMultiply(self.H_w, newXW)
        self.Wind.Ww = Ww[0][0]
        self.x_w = newXW

        return

    def getWind(self):
        """
        Wrapper function to return the wind state from the module
        :return: windState class
        """
        return self.Wind

    def reset(self):
        """
        Wrapper function that resets the wind model code (but does not reset the model chosen for wind. To
        change the model transfer functions you need to use CreateDrydenTranferFns.
        
        :return: none
        """
        self.Wind = States.windState()
        self.x_u= [[0.0]]
        self.x_v = [[0.0], [0.0]]
        self.x_w = [[0.0], [0.0]]
        return
    
    def setWind(self, windState):
        """
        Wrapper function that allows for injecting constant wind and gust values into the class :param windState:
        class from vehicleStates with inertial constant wind and body gusts
        
        :return: none
        """
        self.Wind = windState
        return