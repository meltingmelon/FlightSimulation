import math
import random
from ece163.Modeling import VehicleAerodynamicsModel
from ece163.Utilities import MatrixMath as mm
from ..Containers import Sensors
from ..Constants import VehiclePhysicalConstants as VPC
from ..Constants import VehicleSensorConstants as VSC


class GaussMarkov:
    def __init__(self, dT=VPC.dT, tau=1e6, eta=0.0):
        """
        Function to initialize the GaussMarkov code that generates
        the exponentially correlated noise which is used for the slowly
        varying bias drif of the gyros as well as the GPS position.
        Implements the noise model characterized by a first order Gauss-Markov process:
        dv/dt = -1/tau v + w, where w is a white noise process with N(0,sigma)

        :param dT: time step [s]
        :param tau: correlation time [s]
        :param sigma: standard deviation of white noise process

        :return none:
        """
        self.dT = dT
        self.tau = tau
        self.sigma = eta
        if (self.sigma == None):
            self.w = 0.0
            self.sigma = 0.0
        else:
            self.w = random.gauss(0, self.sigma)
        self.vnoise_prev = 0.0
        return

    def reset(self):
        """
        Wrapper function that resets the GaussMarkov model

        :return none:
        """
        self.dT = VPC.dT
        self.tau = 1e6
        self.eta = 0.0
        self.vnoise_prev = 0.0
        return

    def update(self, vnoise=None):
        """
        Function that updates the Gauss-Markov process, and returns the updated value
        (as well as updating the internal value that holds until the next call of the function)

        :param vnoise: optional parameter to drive the function with a known value. If none, then
        use random.gauss(0,sigma)

        :return v: new noise value (also updated internally)
        """
        dT = self.dT
        tau = self.tau
        w = self.w
        
        if (vnoise == None):
            v = (math.exp(-1 * (dT / tau)) * w) + w
        else:
            v = (math.exp(-1 * (dT / tau)) * self.vnoise_prev) + vnoise
        
        self.vnoise_prev = v
        return v


class GaussMarkovXYZ:
    def __init__(self, dT=VPC.dT, tauX=1e6, etaX=0.0, tauY=None, etaY=None, tauZ=None, etaZ=None):
        """
        Function to aggregate three Gauss-Markov models into a triplet that returns the X, Y, and Z axes
        of the time-varying drift; if (tau,eta) are None, then will default to the same values for each model

        :param dT: time step [s]
        :param tau: correlation time [s]
        :param sigma: standard deviation of white noise process

        :return none:
        """
        
        if (tauY == None):
            tauY = tauX
        if (etaY == None):
            etaY = etaX
        if (tauZ == None):
            tauZ = tauX
        if (etaZ == None):
            etaZ = etaX

        self.vNoiseX = GaussMarkov(dT, tauX, etaX)
        self.vNoiseY = GaussMarkov(dT, tauY, etaY)
        self.vNoiseZ = GaussMarkov(dT, tauZ, etaZ)
        return

    def reset(self):
        """
        Wrapper function that resets the GaussMarkovXYZ models

        :return none:
        """
        
        self.vNoiseX.reset()
        self.vNoiseY.reset()
        self.vNoiseZ.reset()

        return

    def update(self, vXnoise=None, vYnoise=None, vZnoise=None):
        """
        Function that updates the Gauss-Markov processes, and returns the updated values
        (as well as updating the internal values that holds until the next function call)

        :param vXnoise: optional parameter to drive the X function with a known value. If None,
        then use random.gauss(0,sigma)
        :param vYnoise: optional parameter to drive the Y function with a known value, If None,
        then use random.gauss(0,sigma)
        :param vZnoise: optional parameter to drive the Z function with a known value. If None,
        then us random.gauss(0,sigma)

        :return vX: new noise value for X
        :return vY: new noise value for Y
        :return vZ: new noise value for Z
        """

        vX = self.vNoiseX.update(vXnoise)
        vY = self.vNoiseY.update(vYnoise)
        vZ = self.vNoiseZ.update(vZnoise)

        return vX, vY, vZ


class SensorsModel:
    def __init__(self, aeroModel=VehicleAerodynamicsModel.VehicleAerodynamicsModel(), taugyro=VSC.gyro_tau, etagyro=VSC.gyro_eta, tauGPS=VSC.GPS_tau, etaGPSHorizontal=VSC.GPS_etaHorizontal, etaGPSVertical=VSC.GPS_etaVertical, gpsUpdateHz=VSC.GPS_rate):
        """
        Function to initialize the SensorsModel code. Will contain both the true sensors outputs and the 
        noisy sensor outputs using the noise and biases found in the Constants.VehicleSensorConstants file.
        Biases for the sensors are set at instantiation and not updated further during the code run. 
        All sensor outputs are in Engineering Units (it is assumed that the raw ADC counts to engineering unit
        scaling have already been done). Note that if the biases are set to None at instantiation,
        then they will be set to random values using uniform distribution with the appropriate bias scaling from
        the sensors constants. SensorsModel class keeps the Gauss-Markov models for gyro biases and GPS.

        :param dT: time steo [s]
        :param taugyro: Gauss-Markov time constant for gyro [s]
        :param etagyro: Gauss-Markov processs noise standard deviation [rad/s]
        :param tauGPS: Gauss-Markov time constant for GPS [s]
        :oaram etaGPSHorizontal: Gauss-Markov process noise standard deviation [m]
        :param etaGPSVertical: Gauss-Markov process noise standard deviation [m]
        :param gpsUpdateHz: Update rate for GPS measurements [Hz]

        :return none:
        """

        self.sensorsTrue = Sensors.vehicleSensors()
        self.sensorsNoisy = Sensors.vehicleSensors()
        self.sensorBiases = self.initializeBiases()
        self.sensorSigmas = self.initializeSigmas()
        self.VAM = aeroModel
        dT = self.VAM.vehicleDynamics.dT
        self.updateTicks = 0
        self.gpsUpdateTicks = int(1 / (dT * gpsUpdateHz))

        self.gyroGM = GaussMarkovXYZ(dT, taugyro, etagyro)
        self.gpsGM = GaussMarkovXYZ((1 / gpsUpdateHz), tauGPS, etaGPSHorizontal, tauGPS, etaGPSHorizontal, tauGPS, etaGPSVertical)


        return

    def getSensorsNoisy(self):
        """
        Wrapper function to return the noisy sensor values

        :return sensorNoisy: noisy sensor data (Sensors.vehicleSensors class object)
        """
        return self.sensorsNoisy

    def getSensorsTrue(self):
        """
        Wrapper function to return the true sensor values

        :return sensorsTrue: true sensor values (Sensors.vehicleSensors class object)
        """
        return self.sensorsTrue

    def initializeBiases(self, gyroBias=VSC.gyro_bias, accelBias=VSC.accel_bias, magBias=VSC.mag_bias, baroBias=VSC.baro_bias, pitotBias=VSC.pitot_bias):
        """
        Function to generate the biases for each of the sensors. Biases are set with a  uniform random number
        from -1 to 1 that is then multiplied by the sigma_bias. The biases for all sensors is returned as a
        Sensors.vehicleSensors class. Note that GPS is an unbiased sensor (though noisy), thus all the GPS biases
        are set to 0.0

        :param gyroBias: bias scaling for the gyros [rad/s]
        :param accelBias: bias scaling for the accelerometers [m/s^2]
        :param magBias: bias scaling for the magnetometers [nT]
        :param baroBias: bias scaling for the barometer [N/m^2]
        :param pitotBias: bias scaling for the pitot tube [N/m^2]

        :returns sensorBiases: Sensors.vehicleSensors class object
        """

        gyro_x = random.uniform(-1, 1) * gyroBias
        gyro_y = random.uniform(-1, 1) * gyroBias
        gyro_z = random.uniform(-1, 1) * gyroBias

        accel_x = random.uniform(-1, 1) * accelBias
        accel_y = random.uniform(-1, 1) * accelBias
        accel_z = random.uniform(-1, 1) * accelBias

        mag_x = random.uniform(-1, 1) * magBias
        mag_y = random.uniform(-1, 1) * magBias
        mag_z = random.uniform(-1, 1) * magBias

        baro = random.uniform(-1, 1) * baroBias
        pitot = random.uniform(-1, 1) * pitotBias

        sensorBiases = Sensors.vehicleSensors()

        sensorBiases.gyro_x = gyro_x
        sensorBiases.gyro_y = gyro_y
        sensorBiases.gyro_z = gyro_z

        sensorBiases.accel_x = accel_x
        sensorBiases.accel_y = accel_y
        sensorBiases.accel_z = accel_z

        sensorBiases.mag_x = mag_x
        sensorBiases.mag_y = mag_y
        sensorBiases.mag_z = mag_z

        sensorBiases.baro = baro
        sensorBiases.pitot = pitot

        sensorBiases.gps_n = 0.0
        sensorBiases.gps_e = 0.0
        sensorBiases.gps_alt = 0.0
        sensorBiases.gps_sog = 0.0
        sensorBiases.gps_cog = 0.0

        return sensorBiases

    def initializeSigmas(self, gyroSigma=VSC.gyro_sigma, accelSigma=VSC.accel_sigma, magSigma=VSC.mag_sigma, baroSigma=VSC.baro_sigma, pitotSigma=VSC.pitot_sigma, gpsSigmaHorizontal=VSC.GPS_sigmaHorizontal, gpsSigmaVertical=VSC.GPS_sigmaVertical, gpsSigmaSOG=VSC.GPS_sigmaSOG, gpsSigmaCOG=VSC.GPS_sigmaCOG):
        """
        Function to gather all of the white noise standard deviations into a single vehicleSensor class object. These will be used as the input to generating the white
        noise added to each sensor when generating the noisy sensor data

        :param gyroSigma: gyro white noise [rad/s]
        :param accelSigma: accelerometer white noise [m/s^2]
        :param magSigma: magnetometer white noise [nT]
        :param baroSigma: barometer white noise [N/m]
        :param pitotSigma: airspeed white noise [N/m]
        :param gpsSigmaHorizontal: GPS horizontal white noise [m]
        :param gpsSigmaVertical: GPS vertical white noise [m]
        :param gpsSigmaSOG: GPS Speed over ground white noise [m/s]
        :param gpsSigmaCOG: GPS Course over ground white noise, nominal [rad]

        :returns sensorSigmas: Sensors.vehicleSensors class object
        """

        sensorSigmas = Sensors.vehicleSensors()

        sensorSigmas.gyro_x = gyroSigma
        sensorSigmas.gyro_y = gyroSigma
        sensorSigmas.gyro_z = gyroSigma

        sensorSigmas.accel_x = accelSigma
        sensorSigmas.accel_y = accelSigma
        sensorSigmas.accel_z = accelSigma

        sensorSigmas.mag_x = magSigma
        sensorSigmas.mag_y = magSigma
        sensorSigmas.mag_z = magSigma

        sensorSigmas.baro = baroSigma
        sensorSigmas.pitot = pitotSigma

        sensorSigmas.gps_n = gpsSigmaHorizontal
        sensorSigmas.gps_e = gpsSigmaHorizontal
        sensorSigmas.gps_alt = gpsSigmaVertical
        sensorSigmas.gps_sog = gpsSigmaSOG
        sensorSigmas.gps_cog = gpsSigmaCOG

        return sensorSigmas

    def reset(self):
        """
        Function reset the module to run again. Should reset the Gauss-Markov models, re-initialize the sensor biases,
        and reset the sensors true and noisy to prestine conditions.

        :returns none:
        """

        self.sensorsNoisy  = Sensors.vehicleSensors()
        self.sensorsTrue = Sensors.vehicleSensors()
        self.sensorSigmas = self.initializeSigmas()
        self.sensorBiases = self.initializeBiases()

        self.gyroGM.reset()
        self.gpsGM.reset()
        self.updateTicks = 0

        return

    def updateAccelsTrue(self, state, dot):
        """
        Function to update the accelerometer sensor. Will be called within the updateSensors function

        :param state: States.vehicleState class object. Current vehicle state
        :param dot: States.vehicleState class object. Current derivative.

        :returns accel_x,accel_y,accel_z: body frame ref force [m/s^2]
        """


        accel_x = dot.u + (state.q * state.w) - (state.r * state.v) + (VPC.g0 * math.sin(state.pitch))
        accel_y = dot.v + (state.r * state.u) - (state.p * state.w) - (VPC.g0 * math.cos(state.pitch) * math.sin(state.roll))
        accel_z = dot.w + (state.p * state.v) - (state.q * state.u) - (VPC.g0 * math.cos(state.pitch) * math.cos(state.roll))

        return accel_x, accel_y, accel_z

    def updateGPSTrue(self, state, dot):
        """
        Function to update the GPS sensor state (this will be called to update the GPS data from the state
        and the derivative) at the required rate. Note that GPS reports back altitude as + above mean sea level

        :param state: class States.vehicleState, current vehicle state
        :param dot: class States.vehicleState, current state derivative

        :returns gps_n [North - m], gps_e [East - m], gps_alt [Altitude - m], gps_SOG [Speed over ground, m/s],
        gps_COG [course over ground, rad]
        """
        gps_n = state.pn
        gps_e = state.pe
        gps_alt = -1 * state.pd

        gps_SOG = math.hypot(state.u, state.v, state.w)
        gps_COG = math.atan2(dot.pe, dot.pn)

        return gps_n, gps_e, gps_alt, gps_SOG, gps_COG

    def updateGyrosTrue(self, state):
        """
        Function to update the rate gyro sensor. Will be called within the update-Sensors functions.

        :param state: class States.vehicleState

        :return gyro_x,gyro_y,gyro_z: body frame rotation rates [rad/s]
        """

        gyro_x = state.p
        gyro_y = state.q
        gyro_z = state.r

        return gyro_x, gyro_y, gyro_z

    def updateMagsTrue(self, state):
        """
        Function to update the magnetometer sensor. Will be called within in the updateSensors function.

        :param state: class States.vehicleState, current vehicle state

        :returns mag_x, mag_y, mag_z: body frame magnetic field [nT]
        """

        R = state.R
        magField = VSC.magfield

        magHeading = mm.matrixMultiply(R, magField)

        mag_x = magHeading[0][0]
        mag_y = magHeading[1][0]
        mag_z = magHeading[2][0]

        return mag_x, mag_y, mag_z

    def updatePressureSensorsTrue(self, state):
        """
        Function to update the pressure sensors onboard the aircraft.
        Will be called within the updateSensors functions. The two pressure sensors are static pressure (barometer)
        and dynamic pressure (pitot tube). The barometric pressure is referenced off of the ground static pressure
        in VehicleSensorConstants at Pground.

        :param state: class States.vehicleState, current vehicle state

        :return baro, pitot: in [N/m^2]
        """
        rho = VPC.rho
        Va = state.Va
        h = (-1 * state.pd)
        g = VPC.g0

        baro = (-1 * (rho * g * h)) + VSC.Pground

        pitot = (rho * (Va ** 2)) / 2

        return baro, pitot

    def updateSensorsNoisy(self, trueSensors=Sensors.vehicleSensors(), noisySensors=Sensors.vehicleSensors(), sensorBiases=Sensors.vehicleSensors(), sensorSigmas=Sensors.vehicleSensors()):
        """
        Function to generate the noisy sensor data given the true sensor readings, the biases, and the sigmas for the
        white noise on each sensor. The gauss markov models for the gyro biases and GPS positions are updated
        here. The GPS COG white noise must be scaled by the ratio of VPC.initialSpeed / actual ground speed.
        GPS is only updated if the correct number of ticks have gone by to indicate that a new GPS measurement
        should be generated. The GPS COG must be limited to within +/- PI. If no GPS update has occurred, then
        the values for the GPS sensors should be copied from the noisySensors input to the output.

        :param trueSensors: Sensors.vehicleSensors class, true values (no noise)
        :param noisySensors: Sensors.vehicleSensors class, previous noisy sensor values
        :param sensorBiases: Sensors.vehicleSensors class, fixed biases for each sensor
        :param sensorSigmas: Sensors.vehicleSensors class, standard deviation of white noise on each sensor.
        """

        sn = Sensors.vehicleSensors()

        #Sigmas
        sigmaAccelX = sensorSigmas.accel_x
        sigmaAccelY = sensorSigmas.accel_y
        sigmaAccelZ = sensorSigmas.accel_z

        sigmaGyroX = sensorSigmas.gyro_x
        sigmaGyroY = sensorSigmas.gyro_y
        sigmaGyroZ = sensorSigmas.gyro_z

        sigmaMagX = sensorSigmas.mag_x
        sigmaMagY = sensorSigmas.mag_y
        sigmaMagZ = sensorSigmas.mag_z

        sigmaBaro = sensorSigmas.baro
        sigmaPitot = sensorSigmas.pitot

        sigmaGPSN = sensorSigmas.gps_n
        sigmaGPSE = sensorSigmas.gps_e
        sigmaGPSAlt = sensorSigmas.gps_alt
        sigmaSOG = sensorSigmas.gps_sog
        sigmaCOG = sensorSigmas.gps_cog

        #Biases
        biasAccelX = sensorBiases.accel_x
        biasAccelY = sensorBiases.accel_y
        biasAccelZ = sensorBiases.accel_z

        biasGyroX = sensorBiases.gyro_x
        biasGyroY = sensorBiases.gyro_y
        biasGyroZ = sensorBiases.gyro_z

        biasMagX = sensorBiases.mag_x
        biasMagY = sensorBiases.mag_y
        biasMagZ = sensorBiases.mag_z

        biasBaro = sensorBiases.baro
        biasPitot = sensorBiases.pitot

        # Gauss Markov
        gyroGaussX, gyroGaussY, gyroGaussZ = self.gyroGM.update()

        #Noisy Signals
        sn.accel_x = trueSensors.accel_x + biasAccelX + random.gauss(0, sigmaAccelX)
        sn.accel_y = trueSensors.accel_y + biasAccelY + random.gauss(0, sigmaAccelY)
        sn.accel_z = trueSensors.accel_z + biasAccelZ + random.gauss(0, sigmaAccelZ)

        sn.gyro_x = trueSensors.gyro_x + biasGyroX + gyroGaussX + random.gauss(0, sigmaGyroX)
        sn.gyro_y = trueSensors.gyro_y + biasGyroY + gyroGaussY + random.gauss(0, sigmaGyroY)
        sn.gyro_z = trueSensors.gyro_z + biasGyroZ + gyroGaussZ + random.gauss(0, sigmaGyroZ)

        sn.mag_x = trueSensors.mag_x + biasMagX + random.gauss(0, sigmaMagX)
        sn.mag_y = trueSensors.mag_y + biasMagY + random.gauss(0, sigmaMagY)
        sn.mag_z = trueSensors.mag_z + biasMagZ + random.gauss(0, sigmaMagZ)

        sn.baro = trueSensors.baro  + biasBaro + random.gauss(0, sigmaBaro)
        sn.pitot = trueSensors.pitot + biasPitot + random.gauss(0, sigmaPitot)

        if ((self.updateTicks % self.gpsUpdateTicks) == 0):
            gpsBn, gpsBe, gpsAlt = self.gpsGM.update()
            sn.gps_n = trueSensors.gps_n + gpsBn + random.gauss(0, sigmaGPSN)
            sn.gps_e = trueSensors.gps_e + gpsBe + random.gauss(0, sigmaGPSE)
            sn.gps_alt = trueSensors.gps_alt + gpsAlt + random.gauss(0,sigmaGPSAlt)
            sn.gps_sog = trueSensors.gps_sog + random.gauss(0, sigmaSOG)

            if (math.isclose(0, noisySensors.gps_sog)):
                sn.gps_cog = trueSensors.gps_cog + random.gauss(0, (sigmaCOG * 100))
            else:
                sn.gps_cog = trueSensors.gps_cog + random.gauss(0,((sensorSigmas.gps_cog * VPC.InitialSpeed) / trueSensors.gps_sog))
            # wrap cog between +- pi
            if (sn.gps_cog > math.pi):
                sn.gps_cog = math.pi
            if (sn.gps_cog < -1 * math.pi):
                sn.gps_cog = -1 * math.pi

        else:
            #Keep all the same
            sn.gps_n = noisySensors.gps_n
            sn.gps_e = noisySensors.gps_e
            sn.gps_alt = noisySensors.gps_alt
            sn.gps_sog = noisySensors.gps_sog
            sn.gps_cog = noisySensors.gps_cog

        return sn

    def updateSensorsTrue(self, preTrueSensors, state, dot):
        """
        Function to generate the true sensors given the current state and
        state derivative. Sensor suite is 3-axis accelerometer, 3-axis rate gyros, 3-axis magnetometers, a barometric
        altimeter, a pitot airspeed, and GPS with an update rate specified in the VehicleSensorConstants file. For
        the GPS update, the previous value is returned until a new update occurs.

        :param state: class States.vehicleState, current vehicle state
        :param dot: class States.vehicleState, current state derivative

        :returns sensorsTrue: Sensors.vehicleSensors class
        """

        sensorsTrue=Sensors.vehicleSensors()

        sensorsTrue.accel_x, sensorsTrue.accel_y, sensorsTrue.accel_z = self.updateAccelsTrue(state, dot)
        sensorsTrue.gyro_x, sensorsTrue.gyro_y, sensorsTrue.gyro_z = self.updateGyrosTrue(state)
        sensorsTrue.mag_x, sensorsTrue.mag_y, sensorsTrue.mag_z = self.updateMagsTrue(state)


        sensorsTrue.baro, sensorsTrue.pitot = self.updatePressureSensorsTrue(state)

        if ((self.updateTicks % self.gpsUpdateTicks) == 0):
            sensorsTrue.gps_n, sensorsTrue.gps_e, sensorsTrue.gps_alt,  sensorsTrue.gps_sog, sensorsTrue.gps_cog = self.updateGPSTrue(state, dot)
        else:
            sensorsTrue.gps_n = preTrueSensors.gps_n
            sensorsTrue.gps_e = preTrueSensors.gps_e
            sensorsTrue.gps_alt = preTrueSensors.gps_alt
            sensorsTrue.gps_sog = preTrueSensors.gps_sog
            sensorsTrue.gps_cog = preTrueSensors.gps_cog

        return sensorsTrue

    def update(self):
        state = self.VAM.vehicleDynamics.state
        dot = self.VAM.vehicleDynamics.dot
        prevSensorTrue = self.sensorsTrue
        self.sensorsTrue = self.updateSensorsTrue(prevSensorTrue, state, dot)
        self.sensorsNoisy = self.updateSensorsNoisy(self.sensorsTrue, self.sensorsNoisy, self.sensorBiases, self.sensorSigmas)

        self.updateTicks  += 1

        return
