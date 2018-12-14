class Node:

    def __init__(self, id, x, y, temp):
        self.id = id
        self.x = x
        self.y = y
        self.temp = temp

    def print(self):
        print("ID: {} ({}; {}) temp: {}".format(self.id, self.x, self.y, self.temp))


class Element:

    def __init__(self, id, node1, node2, node3, node4, k):
        self.id = id
        self.nodes = [node1, node2, node3, node4]
        self.k = k

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
