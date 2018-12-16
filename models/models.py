from numpy import  *


class Node:

    def __init__(self, id, x, y, temp):
        self.id = id
        self.x = x
        self.y = y
        self.temp = temp

    def print(self):
        print("ID: {} ({}; {}) temp: {}".format(self.id, self.x, self.y, self.temp))


class Element:

    def __init__(self, id, node1, node2, node3, node4, surface_indexes, k):
        self.id = id
        self.nodes = [node1, node2, node3, node4]
        self.heated_surfaces_indexes = surface_indexes
        self.k = k
        self.H_matrix = zeros([4, 4])
        self.C_matrix = zeros([4, 4])
        self.P_vector = zeros(4)

    def print(self):
        print("ID: {}; Nodes[{}, {}, {}, {}]; k = {}"
              .format(self.id, self.nodes[0].id, self.nodes[1].id,
                      self.nodes[2].id, self.nodes[3].id, self.k))


class Grid:

    def __init__(self):
        self.nodes = []
        self.elements = []

    def __del__(self):
        for node in self.nodes:
            del node
        for element in self.elements:
            del element

    def print(self):
        for node in self.nodes:
            node.print()

        print()

        for element in self.elements:
            element.print()


class PointKsiEta:
    def __init__(self, ksi, eta):
        self.ksi = ksi
        self.eta = eta

