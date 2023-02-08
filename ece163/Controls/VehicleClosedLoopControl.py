import math
import sys
import pickle
import enum
import ece163.Containers.Inputs as Inputs
import ece163.Containers.States as States
import ece163.Containers.Controls as Controls
import ece163.Constants.VehiclePhysicalConstants as VPC
import ece163.Modeling.VehicleAerodynamicsModel as VehicleAerodynamicsModule


class PDControl:
    def __init__(self, kp=0.0, kd=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        """
        Function which implement the PD control with saturation where the derivative is available as a separate input to the function.
        The output is: u = u_ref + Kp * error - Kd * dot{error} limited between lowLimit and highLimit.

        :param kp: proportional gain
        :param kd: derivative gain
        :param trim: trim output (added the loop computed output)
        :param lowLimit: lower limit to saturate control
        :param highLimit: upper limit to saturate control

        :return: none
        """
        self.kp = kp
        self.kd = kd
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit

        return

    def Update(self, command=0.0, current=0.0, derivative=0.0):
        """
        Calculates the output of the PD loop given the gains and limits from instantiation,
        and using the command, actual, and derivative inputs. Output is limited to between
        lowLimit and highLimit from instantiation.

        :param command: reference command
        :param current: actual output (or sensor)
        :param derivative: derivative of the output sensor

        :return u: control limited to saturation bounds.
        """

        e = command - current
        u = (self.kp * e) - (self.kd * derivative) + self.trim

        if (u > self.highLimit):
            u = self.highLimit
        elif (u < self.lowLimit):
            u = self.lowLimit

        return u

    def setPDGains(self, kp=0.0, kd=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        """
        Function to set the gains for the PD control block(including the trim output and the limit)

        :param kp: proportional gain
        :param kd: derivative gain
        :param trim: trim output (Added the loop computed output)
        :param lowLimit: lower limit of saturate control
        :param highLimit: upper limit of saturate control

        :return: none
        """

        self.kp = kp
        self.kd = kd
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit

        return


class PIControl:
    def __init__(self, dT=0.01, kp=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        """
        Function which implement the PI control with saturation where the integrator has both a reset
        and an anti-windup such that when output saturates, the integration is undone and the output
        forced the output to the limit. The output is: u = u_ref + Kp * error + Ki*integral{error} limited
        between lowLimit and highLimit

        :param dT: time step[s], required for integration
        :param kp: proportional gain
        :param ki: integral gain
        :param trim: trim input
        :param lowLimit: low saturation limit
        :param highLimit: high saturation limit

        :return: none
        """
        self.dT = dT
        self.kp = kp
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        self.accumulator = 0.0
        self.prevError = 0.0

        return

    def Update(self, command=0.0, current=0.0):
        """
        Calculates the output of the PI loop given the gains and limits from instantiation, and using the command and current or actual inputs.
        Output is limited to between lowLimit and highLimit from instantiation. Integration for the integral state is done using trapezoidal
        integration, and anti-windup is implemented such that if the output is out of limits, the integral state is not updated (no additional error
        accumulation).

        :param command: reference command
        :param current: current output or sensor

        :return u: output limited to saturation bounds
        """
        e = command - current
        self.accumulator += (1 / 2) * (self.dT * (e + self.prevError))
        u = (self.kp * e) + (self.ki * self.accumulator) + self.trim

        if (u > self.highLimit):
            u = self.highLimit
            self.accumulator -= (1 / 2) * self.dT * (e + self.prevError)
        elif (u < self.lowLimit):
            u = self.lowLimit
            self.accumulator -= (1 / 2) * self.dT * (e + self.prevError)

        self.prevError = e
        return u

    def resetIntegrator(self):
        """
        Function to reset the integration state to zero, used when switching modes or otherwise resetting the integral state.

        :return: none
        """
        self.accumulator = 0.0
        self.prevError = 0.0
        return

    def setPIGains(self, dT=VPC.dT, kp=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        """
        Function to set the gains for the PI control block (including the trim output and the limits)

        :param dT: time step [s], required for integration
        :param kp: proportional gain
        :param ki: integral gain
        :param trim: trim input
        :param lowLimit: low saturation limit
        :param highLimit: high saturation limit

        :return: none
        """

        self.dT = dT
        self.kp = kp
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit

        return


class PIDControl:
    def __init__(self, dT=VPC.dT, kp=0.0, kd=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        """
         Functions which implement the PID control with saturation where the integrator has both a reset and an anti windup sud
         what when output saturates, the integration is undone and the output forced the output to the limit. 
         Function assumes that physical derivative is available that physical derivative is available (e.g. roll and p), not a numerically derived one.

         :param dT: time step [s], required for integration
         :param kp: proportional gain
         :param kd: derivative gain
         :param ki: integral gain
         :param trim: trim input
         :param lowLimit: low saturation limit
         :param highLimit: high saturation limit

         :return: none
        """

        self.dT = dT
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit

        self.accumulator = 0.0
        self.prevError = 0.0

        return

    def Update(self, command=0.0, current=0.0, derivative=0.0):
        """
        Calculates the output of the PID loop given the gains and limits from instantiation, and using the command and current or actual inputs. Output is
        limited to between lowLimit and highLimit from instantiation. Integration for the integral state is done
        using trapezoidal integration, and anti-windup is implemented such that if the output is out of limits, the
        integral state is not updated (no additional error accumulation).

        :param command: reference command
        :param current: current error output or sensor
        :param derivative: derivative of the output or sensor

        :return u: output limited to saturation bounds.
        """
        e = command - current

        self.accumulator += (1 / 2) * self.dT * (e + self.prevError)
        u = (self.kp * e) + (-1 * self.kd * derivative) + (self.ki * self.accumulator) + self.trim

        if (u > self.highLimit):
            u = self.highLimit
            self.accumulator -= (1 / 2) * self.dT * (e + self.prevError)
        elif (u < self.lowLimit):
            u = self.lowLimit
            self.accumulator -= (1 / 2) * self.dT * (e + self.prevError)

        self.prevError = e
        return u

    def resetIntegrator(self):
        """
        Function to reset the integration state to zero, used when switching modes or otherwise resetting the integral state.
        """
        self.accumulator = 0.0
        self.prevError = 0.0

        return

    def setPIDGains(self, dT=VPC.dT, kp=0.0, kd=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        """
        Function to set the gains for the PI control block(including the trim output and the limits)

        :param dT: time step [s], required for integration
        :param kp: proportional gain
        :param kd: derivative gain
        :param ki: integral gain
        :param trim: trim input
        :param lowLimit: low saturation limit
        :param highLimit: high saturation limit

        :return: none
        """

        self.dT = dT
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit

        return


class VehicleClosedLoopControl:
    def __init__(self, dT=VPC.dT):
        """
        Contains the required state for the altitude hold state machine from the enumeration class in Controls.AltitudeStates

        :param dT: time step [s]

        :return: none
        """

        self.controlGains = Controls.controlGains()
        self.trimInputs = Inputs.controlInputs()
        self.VehicleControlSurfaces = Inputs.controlInputs()

        self.VAM = VehicleAerodynamicsModule.VehicleAerodynamicsModel()
        self.VAM.vehicleDynamics.dT = dT
        self.dT = self.VAM.vehicleDynamics.dT

        self.aileronFromRoll = PIDControl()
        self.rollFromCourse = PIControl()
        self.rudderFromSideslip = PIControl()
        self.elevatorFromPitch = PDControl()
        self.throttleFromAirspeed = PIControl()
        self.pitchFromAirspeed = PIControl()
        self.pitchFromAltitude = PIControl()
        self.climbMode = Controls.AltitudeStates.HOLDING

        return

    def Update(self, referenceCommands=Controls.referenceCommands()):
        """
        Function that implements the full closed loop controls using the command inputs of airspeed, altitude, and course.
        Calls all of the submodules below it to implement the full flight dynamics under closed loop control.
        Note that the internal commandPitch and commandedRoll of the reference commands are altered by this function, and this
        is used to track the reference commands by the simulator.

        You will need to add or subtract 360 degree (2 * pi) from you internal course when the error is outside
        of +/- 180 degrees(pi radians). There is never a course error of more than 180 degrees, so if you do see
        such an error, it is because signed angles and reciprocal headings(e.g.: -180 and +180 are pointed in the same direction).

         :param referenceCommands: high level autopilot commands (note: altered by function)
         :return: none 
         pseudocode:
        courseError = ref.commandedCourse - VAM,VD.state.chi
        if courseError >= pi or 

        rollCommand = rollFromCourse.Update(commandedCourse, currentCourse)
        aileronCommand = aileronFromRoll.Update(roll command, roll from state, p)
        rudder = rudderFromSideslip = Update(0.0, beta from state)

        altitude = -pd
        if altitude > (ref.commandedAltitude + VPC.altitudeHoldZone):
            if self.climbMode is not control.AltState.descending:
            self.climbMode = descending
            pitchFromAirspeed.resetIntegrator()alright so I think I'd like to see if 
        commandedThrottle = VPC.minControlThrottle
        commandedPitch = pitchFromAirspeed.Update(ref.commandedAirspeed, airspeed from state)

        elif altitude < (ref.commandedAltitude - VPC.altitudeHoldZone):
            do the opposite

        else:
            if climbMode != holding
                climbMode = holding
            throttleCommand = throttleFromAirspeed.update(ref.commandAirspeed, Va from state)
            pitchCommand = pitchFromAltitude.Update(ref cmdAlt, altitude)

        elevatorCommand = self.elevatorfromPitch.Update(pitchCommand, pitch from state)

        self.ControlSurfaceOutputs.Throttle
        .elevator
        .aileron
        .rudder
        (set to all of these variables)
        refcommandedPitch = pitchCommand 
        refroll = rollCommand
        """
        courseError = referenceCommands.commandedCourse - self.VAM.vehicleDynamics.state.chi

        if (courseError >= math.pi):
            self.VAM.vehicleDynamics.state.chi += (2 * math.pi)

        if (courseError <= (-1 * math.pi)):
            self.VAM.vehicleDynamics.state.chi -= (2 * math.pi)

        ########### Lateral/Directional #############
        courseCommand = referenceCommands.commandedCourse
        courseCurrent = self.VAM.vehicleDynamics.state.chi

        # Aileron parameters
        rollCommand = self.rollFromCourse.Update(courseCommand, courseCurrent)
        rollCurrent = self.VAM.vehicleDynamics.state.roll
        p = self.VAM.vehicleDynamics.state.p

        aileronCommand = self.aileronFromRoll.Update(rollCommand, rollCurrent, p)

        # Rudder parameters
        betaCommand = 0.0
        betaCurrent = self.VAM.vehicleDynamics.state.beta

        # Altitude parameters
        altitudeCommand = referenceCommands.commandedAltitude
        altitudeCurrent = -1 * self.VAM.vehicleDynamics.state.pd

        # Airspeed parameters
        airspeedCommand = referenceCommands.commandedAirspeed
        airspeedCurrent = self.VAM.vehicleDynamics.state.Va
        rudderCommand = self.rudderFromSideslip.Update(0.0, betaCurrent)


        if (altitudeCurrent > (altitudeCommand + VPC.altitudeHoldZone)):
            if (self.climbMode != Controls.AltitudeStates.DESCENDING):
                self.climbMode = Controls.AltitudeStates.DESCENDING
                self.pitchFromAirspeed.resetIntegrator()
                self.pitchFromAltitude.resetIntegrator()
            throttleCommand = VPC.minControls.Throttle
            pitchCommand = self.pitchFromAirspeed.Update(airspeedCommand, airspeedCurrent)

        elif (altitudeCurrent < (altitudeCommand - VPC.altitudeHoldZone)):
            if(self.climbMode != Controls.AltitudeStates.CLIMBING):
                self.climbMode = Controls.AltitudeStates.CLIMBING
                self.pitchFromAirspeed.resetIntegrator()
                self.pitchFromAltitude.resetIntegrator()
            throttleCommand = VPC.maxControls.Throttle
            pitchCommand = self.pitchFromAirspeed.Update(airspeedCommand, airspeedCurrent)

        else:
            if (self.climbMode != Controls.AltitudeStates.HOLDING):
                self.climbMode = Controls.AltitudeStates.HOLDING
                self.pitchFromAirspeed.resetIntegrator()
                self.pitchFromAltitude.resetIntegrator()
            throttleCommand = self.throttleFromAirspeed.Update(airspeedCommand, airspeedCurrent)
            pitchCommand = self.pitchFromAltitude.Update(altitudeCommand, altitudeCurrent)

        pitchCurrent = self.VAM.vehicleDynamics.state.pitch
        q = self.VAM.vehicleDynamics.state.q
        elevatorCommand = self.elevatorFromPitch.Update(pitchCommand, pitchCurrent, q)

        self.VehicleControlSurfaces.Aileron = aileronCommand
        self.VehicleControlSurfaces.Elevator = elevatorCommand
        self.VehicleControlSurfaces.Rudder = rudderCommand
        self.VehicleControlSurfaces.Throttle = throttleCommand

        referenceCommands.commandedPitch = pitchCommand
        referenceCommands.commandedRoll = rollCommand

        self.VAM.Update(self.VehicleControlSurfaces)
        return

    def getControlGains(self):
        """
        Wrapper function to extract control gains from the class.

        :return controlGains: 
        """
        return self.controlGains

    def getVehicleAerodynamicsModel(self):
        """
        Wrapper function to extract the internal VehicleAerodynamicsModel in order to access the various functions
        that are associated with the Aero model (such as setting and getting the wind state and model)

        :return VehicleAerodynamicsModel:
        """

        return self.VAM

    def getVehicleControlSurfaces(self):
        """
        Wrapper function to extract control outputs(Throttle, Aileron, Elevator, Rudder) from the class.

        :return controlInputs:
        """

        return self.VehicleControlSurfaces

    def getVehicleState(self):
        """
        Wrapper function to extract vehicle state from the class

        :return vehicleState:
        """

        return self.VAM.vehicleDynamics.state

    def reset(self):
        """
        Resets the module to run again. Does not overwrite control gains, but does reset the integral states of all
        of the PI control loops.

        :return: none
        """
        self.VAM.reset()
        self.aileronFromRoll.resetIntegrator()
        self.rollFromCourse.resetIntegrator()
        self.rudderFromSideslip.resetIntegrator()
        self.throttleFromAirspeed.resetIntegrator()
        self.pitchFromAirspeed.resetIntegrator()
        self.pitchFromAltitude.resetIntegrator()

        return

    def getTrimInputs(self):
        return self.trimInputs

    def setControlGains(self, controlGains=Controls.controlGains()):
        """
        Function to set all the gains from the controlGains previously computed to the correct places within the various
        control loops to do the successive looop closure (see Beard Chapter 6). Control loop limits are taken from the
        VehiclePhysicalConstants file, trim inputs are taken from self.trimInputs.

        :param controlGains: controlGains class, has the PID loop gains for each loop

        :return: none
        """
        self.controlGains = controlGains

        courseMax = math.radians(VPC.bankAngleLimit)
        courseMin = -1 * math.radians(VPC.bankAngleLimit)

        pitchMax = math.radians(VPC.pitchAngleLimit)
        pitchMin = -1 * math.radians(VPC.pitchAngleLimit)

        aileronMin = VPC.minControls.Aileron
        aileronMax = VPC.maxControls.Aileron

        elevatorMin = VPC.minControls.Elevator
        elevatorMax = VPC.maxControls.Elevator

        rudderMin = VPC.minControls.Rudder
        rudderMax = VPC.maxControls.Rudder

        throttleMin = VPC.minControls.Throttle
        throttleMax = VPC.maxControls.Throttle

        dT = self.dT

        kp_roll = self.controlGains.kp_roll
        kd_roll = self.controlGains.kd_roll
        ki_roll = self.controlGains.ki_roll

        kp_sideslip = self.controlGains.kp_sideslip
        ki_sideslip = self.controlGains.ki_sideslip

        kp_course = self.controlGains.kp_course
        ki_course = self.controlGains.ki_course

        kp_pitch = self.controlGains.kp_pitch
        kd_pitch = self.controlGains.kd_pitch

        kp_altitude = self.controlGains.kp_altitude
        ki_altitude = self.controlGains.ki_altitude

        kp_SpeedfromThrottle = self.controlGains.kp_SpeedfromThrottle
        ki_SpeedfromThrottle = self.controlGains.ki_SpeedfromThrottle

        kp_SpeedfromElevator = self.controlGains.kp_SpeedfromElevator
        ki_SpeedfromElevator = self.controlGains.ki_SpeedfromElevator

        trimInputs = self.trimInputs

        self.aileronFromRoll.setPIDGains(dT, kp_roll, kd_roll, ki_roll, trimInputs.Aileron, aileronMin, aileronMax)
        self.rollFromCourse.setPIGains(dT, kp_course, ki_course, 0.0, courseMin, courseMax)
        self.rudderFromSideslip.setPIGains(dT, kp_sideslip, ki_sideslip, trimInputs.Rudder, rudderMin, rudderMax)
        self.elevatorFromPitch.setPDGains(kp_pitch, kd_pitch, trimInputs.Elevator, elevatorMin, elevatorMax)
        self.throttleFromAirspeed.setPIGains(dT, kp_SpeedfromThrottle, ki_SpeedfromThrottle, trimInputs.Throttle, throttleMin, throttleMax)
        self.pitchFromAirspeed.setPIGains(dT, kp_SpeedfromElevator, ki_SpeedfromElevator, 0.0,  pitchMin, pitchMax)
        self.pitchFromAltitude.setPIGains(dT, kp_altitude, ki_altitude, 0.0, pitchMin, pitchMax)

        return

    def setTrimInputs(self, trimInputs=Inputs.controlInputs(Throttle=0.5, Aileron=0.0, Elevator=0.0, Rudder=0.0)):
        """
        Wrapper function to inject the trim inputs into the class

        :param trimInputs: from Inputs.controlInputs

        :return: none
        """

        self.trimInputs = trimInputs

        return

    def setVehicleState(self, state):
        """
        Wrapper function to inject vehicle state into the class.

        :param state: vehicleState class

        :return: none
        """

        self.VAM.vehicleDynamics.state = state

        return
