import numpy as np
import ece163.Utilities.MatrixMath as mm
import random

class Waypoints:
    def __init__(self, waypoint_type = 'straight line'):
        """
        Waypoints for Path Planning

        :param waypoint_changed: flag to indicate waypoints recently changed (set by planner)
        :param waypoint_manager_request: flag to indicate the waypoint manager needs new points (we may not use this)
        :param type: type of waypoint following. (straight line, fillets (curves), dubins path)
        :para num_waypoints: current number of valid waypoints
        :param ned: ned coordinates of waypoints
        :param airspeed: airspeed command along the waypoints
        :param course: desired course at each waypoint (Dubins only)
        :param cost: cost at each node
        :param parent: index of the parent node
        :param connect_goal: node connect to the goal

        :return none:
        """
        self.waypoint_changed = True
        self.waypoint_manager_request = True

        self.type = waypoint_type
        self.num_waypoints = 0

        self.ned = np.array([[], [], []])
        self.airspeed = np.array([])
        self.course = np.array([])
        self.cost = np.array([])
        sel.parent = np.array([])
        self.connect_goal = np.array([])

        #self.ned = []
        #self.airspeed = []
        #self.course = []
        #self.parent = []
        #self.connect_goal = []

        return

    def add(self, ned, airspeed, course, cost, parent, connect_goal):
        """
        Adds new values into each waypoint array.

        :param ned: new ned coords
        :param airspeed: new airspeed
        :param course: new course
        :param cost: new cost value
        :param parent: new parent value
        :param connect_goal: new connect_goal

        :return none:
        """
        self.num_waypoints = self.num_waypoints + 1
        self.ned = np.append(self.ned,ned,axis=1)
        self.airspeed = np.append(self,airspeed, airspeed)
        self.course = np.append(self.course, course)
        self.cost = np.append(self.cost,cost)
        self.parent = np.append(self.parent,parent)
        self.connect_goal = np.append(self.connect_goal, connect_goal)
        
        return
        
