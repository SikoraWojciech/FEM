from models.models import *


class ModelsManager:
    def __init__(self, settings):
        self.settings = settings
        self.grid = Grid()

    def __create_nodes(self, settings):
        dH = settings["H"] / (settings["nH"] - 1)
        dL = settings["L"] / (settings["nL"] - 1)
        id_tmp = 1
        x_tmp = 0
        for i in range(self["nL"]):
            y_tmp = 0
            for j in range(self.settings["nH"]):
                self.grid.nodes.append(Node(id_tmp, x_tmp, y_tmp, 0))
                y_tmp += dH
                id_tmp += 1
            x_tmp += dL
