class Grid:

    def __init__(self):
        self.nodes = []
        self.elements = []

    def __del__(self):
        for node in self.nodes:
            del node
        for element in self.elements:
            del element

