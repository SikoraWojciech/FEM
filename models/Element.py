class Element:

    def __init__(self, id, node1, node2, node3, node4, k):
        self.id = id
        self.nodes = [node1, node2, node3, node4]
        self.k = k  # wspolczynnik przewodzenia

    def print_element(self):
        print("ID: {}; Nodes[{}, {}, {}, {}]; k = {}"
              .format(self.id, self.nodes[0].id, self.nodes[1].id,
                      self.nodes[2].id, self.nodes[3].id, self.k))
