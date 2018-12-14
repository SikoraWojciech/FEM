class Node:

    def __init__(self, id, x, y, temp):
        self.id = id
        self.x = x
        self.y = y
        self.temp = temp

    def print_node(self):
        print("ID: {} ({}; {}) temp: {}".format(self.id, self.x, self.y, self.temp))

