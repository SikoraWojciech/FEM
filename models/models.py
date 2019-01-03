from numpy import *


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

    def __init__(self, settings):
        self.nodes = []
        self.elements = []
        nodes_count = settings["nH"] * settings["nL"]
        self.H_matrix = zeros((nodes_count, nodes_count))
        self.C_matrix = zeros((nodes_count, nodes_count))
        self.P_vector = zeros(nodes_count)

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

    def agregate_matrices(self):
        for element in self.elements:
            for i_local in range(4):
                for j_local in range(4):
                    i_global = element.nodes[i_local].id - 1
                    j_global = element.nodes[j_local].id - 1

                    val_H_local = element.H_matrix.item((i_local, j_local))
                    val_H_global = self.H_matrix.item((i_global, j_global)) + val_H_local
                    self.H_matrix.itemset((i_global, j_global), val_H_global)

                    val_C_local = element.C_matrix.item((i_local, j_local))
                    val_C_global = self.C_matrix.item((i_global, j_global)) + val_C_local
                    self.C_matrix.itemset((i_global, j_global), val_C_global)

    def agregate_vector(self):
        for element in self.elements:
            for i_local in range(4):
                i_global = element.nodes[i_local].id - 1

                val_P_local = element.P_vector.item(i_local)
                val_P_global = self.P_vector.item(i_global) + val_P_local
                self.P_vector.itemset(i_global, val_P_global)

    def heat(self, settings):
        nodes_count = len(self.nodes)
        dT = settings["time_step"]

        left_part_of_equation = self.H_matrix + (self.C_matrix / dT)

        for interval in range(0, settings["simulation_time"], dT):

            # Wektor t0 - temperatury kazdego z wezlow
            t0_vector = zeros(nodes_count)
            for i in range(nodes_count):
                t0_vector.itemset(i, self.nodes[i].temp)

            right_part_of_equation = -self.P_vector + (self.C_matrix / dT) @ t0_vector

            t1_vector = linalg.solve(left_part_of_equation, right_part_of_equation)

            for i in range(len(t1_vector)):
                self.nodes[i].temp = t1_vector[i]

            print('Interval {}'.format(interval + dT))
            print('=========================================')
            print('P vector ([C]/dT)*{t0} - {P}')
            print('-----------------------------------------')
            print(right_part_of_equation)
            print('-----------------------------------------')
            print('Max: {}'.format(max(t1_vector)))
            print('Min: {}'.format(min(t1_vector)))
            print()


class PointKsiEta:
    def __init__(self, ksi, eta):
        self.ksi = ksi
        self.eta = eta
