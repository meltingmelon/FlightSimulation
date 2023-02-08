import math
from . import MatrixMath

def dcm2Euler(DCM):
    """
    Extracts the Euler angles from the Direction Cosine Matrix

    :param DCM: Rotation matrix [3x3]
    :return: yaw, pitch, roll (radians)
    """
    # Conditional: Clamp DCM[0][2] to +- 1 
    # if value is greater than or less than +-1
    # Implemented based on suggestions from the attitude cheat sheet
    if DCM[0][2] > 1 : DCM[0][2] = 1
    if DCM[0][2] < -1: DCM[0][2] = -1


    # Extract yaw, pitch, roll from DCM using atan2
    # Conversion from Matlab code from Max Dunne in 167,  it is also referenced in the eqn. 3.15, 3.16, and 3.17 in Modeling and Simulation of
    #  Aerospace Vehicles 
    # We can see that these values are just inverses of particular values in the DCM matrix
    yaw = math.atan2(DCM[0][1], DCM[0][0])
    pitch = math.asin(-1*DCM[0][2])
    roll = math.atan2(DCM[1][2], DCM[2][2])

    return yaw, pitch, roll

def euler2DCM(yaw,pitch,roll):
    """
    Creates the Direction Cosine Matrix (DCM) from the Euler angles provided.

    :param yaw: scalar (rotation about inertial down)
    :param pitch: scalar (rotation about intermediate y-axis)
    :param roll: scalar (rotation about body x-axis)
    :return: DCM Matrix
    """

    # Yaw
    DCM_z = [
        [math.cos(yaw), math.sin(yaw), 0],
        [-math.sin(yaw), math.cos(yaw), 0],
        [0, 0, 1]
        ]
    
    # Pitch
    DCM_y = [
        [math.cos(pitch), 0, -math.sin(pitch)],
        [0, 1, 0],
        [math.sin(pitch), 0, math.cos(pitch)]
        ]
    
    # Roll
    DCM_x = [
        [1, 0 , 0], 
        [0, math.cos(roll), math.sin(roll)],
        [0, -math.sin(roll), math.cos(roll)]
        ]

    # 3-2-1, based on equation from Small Unmanned Aircraft eqn 2.5
    DCM = MatrixMath.matrixMultiply(MatrixMath.matrixMultiply(DCM_x,DCM_y),DCM_z)
    return DCM

def ned2enu(points):
    """
    Changes coordinates from North-East-Down (NED) to East-North-Up (ENU).

    :param points: matrix of points in NED [nx3]
    :return: ENU_points same set of [nx3] in ENU coordinates
    """

    # I learned this trick from watching you
    # Creates a new matrix which swaps the column values of NED to ENU
    ENU_points = [[pts[1],pts[0], -pts[2]] for pts in points]

    return ENU_points