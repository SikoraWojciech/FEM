from models.models import *
from shape_functions.matrices import H_matrix_local, C_matrix, P_vector


class ModelsManager:
    def __init__(self, settings):
        self.settings = settings
        self.grid = Grid()

    def __create_nodes(self):
        dH = self.settings["H"] / (self.settings["nH"] - 1)
        dL = self.settings["L"] / (self.settings["nL"] - 1)
        id_tmp = 1
        x_tmp = 0
        for i in range(self.settings["nL"]):
            y_tmp = 0
            for j in range(self.settings["nH"]):
                self.grid.nodes.append(Node(id_tmp, x_tmp, y_tmp, self.settings["initial_temp"]))
                y_tmp += dH
                id_tmp += 1
            x_tmp += dL

    def __create_elements(self):
        nH = self.settings["nH"]
        nL = self.settings["nL"]
        id_tmp = 1
        id_n1_tmp = 0
        for i in range(nL - 1):
            id_n1_tmp = id_n1_tmp + 1
            id_n2_tmp = id_n1_tmp + nH
            id_n3_tmp = id_n2_tmp + 1
            id_n4_tmp = id_n1_tmp + 1

            for j in range(nH - 1):
                element = Element(id_tmp,
                                  self.grid.nodes[id_n1_tmp - 1],
                                  self.grid.nodes[id_n2_tmp - 1],
                                  self.grid.nodes[id_n3_tmp - 1],
                                  self.grid.nodes[id_n4_tmp - 1],
                                  [],
                                  self.settings["k"])
                self.grid.elements.append(element)

                # warunki brzegowe
                if self.grid.nodes[id_n1_tmp - 1].y == 0 \
                        and self.grid.nodes[id_n2_tmp - 1].y == 0:
                    element.heated_surfaces_indexes.append(1)

                if self.grid.nodes[id_n2_tmp - 1].x == self.settings["L"] \
                        and self.grid.nodes[id_n3_tmp - 1].x == self.settings["L"]:
                    element.heated_surfaces_indexes.append(2)

                if self.grid.nodes[id_n3_tmp - 1].y == self.settings["H"] \
                        and self.grid.nodes[id_n4_tmp - 1].y == self.settings["H"]:
                    element.heated_surfaces_indexes.append(3)

                if self.grid.nodes[id_n4_tmp - 1].x == 0 \
                        and self.grid.nodes[id_n1_tmp - 1].x == 0:
                    element.heated_surfaces_indexes.append(4)

                id_tmp += 1
                id_n1_tmp += 1
                id_n2_tmp += 1
                id_n3_tmp += 1
                id_n4_tmp += 1

    def create_grid(self):
        self.__create_nodes()
        self.__create_elements()
        for element in self.grid.elements:
            element.H_matrix = H_matrix_local(element, self.settings["alfa"])
            element.C_matrix = C_matrix(element, self.settings["c"], self.settings["ro"])
            element.P_vector = P_vector(element, self.settings["ambient_temp"], self.settings["alfa"])
