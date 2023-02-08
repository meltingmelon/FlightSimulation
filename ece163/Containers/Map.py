"""
Map Object class

Used to m
"""
import numpy as np
import math
import random
import ece163.Utilities.MatrixMath as mm

class terrainMap:
    def __init__(self, width=2000, height=300, blocks=5, streetWidth=.8):
        """
        Generates random terrain map

        :param map_changed_flag: flag to indicate if map has changed
        :param city_width: the city size [width x width]
        :param num_city_blocks: number of blocks in the city
        :param street_width: the percentage of block that is street
        :param max_building_height: max height of buildings
        :param building_height: array of building heights
        :param building_width: width of buildings (fixed)
        :param building_north: north coordinate of center of buildings
        :param building_east: east coordinate of center of buildings
        """
        self.map_changed_flag = 0
        self.city_width = width
        self.num_city_blocks = blocks
        self.street_width = width / blocks * streetWidth
        self.max_building_height = height
        self.building_height = height * np.random.rand(blocks, blocks)
        # MatrixMath implementation
        #build_array = [[random.uniform(0,1) for i in range(blocks)] for j in range(blocks)]
        #self.building_height = mm.matrixScalarMultiply(height,building_array)
        self.building_width = width / blocks * (1 - streetWidth)

        self.building_north = np.zeros((1,blocks))
        # self.building_north = [[0]] * blocks
        for i in range(blocks):
            self.building_north[0,i] = 0.5 * (width / blocks) * (2 * i + 1)

        self.building_east = np.copy(self.building_north)

        return
        
        
