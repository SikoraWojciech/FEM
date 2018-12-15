from models.models import *


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
                self.grid.nodes.append(Node(id_tmp, x_tmp, y_tmp, 0))
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
                self.grid.elements.append(Element(id_tmp,
                                                  self.grid.nodes[id_n1_tmp - 1],
                                                  self.grid.nodes[id_n2_tmp - 1],
                                                  self.grid.nodes[id_n3_tmp - 1],
                                                  self.grid.nodes[id_n4_tmp - 1],
                                                  self.settings["k"]))
                id_tmp += 1
                id_n1_tmp += 1
                id_n2_tmp += 1
                id_n3_tmp += 1
                id_n4_tmp += 1

    def create_grid(self):
        self.__create_nodes()
        self.__create_elements()
