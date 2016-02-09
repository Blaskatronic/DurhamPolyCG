class NDSparseMatrix:
    def __init__(self):
        self.elements = {}

    def addValue(self, tuple, value):
        self.elements[tuple] = value

    def readValue(self.tuple):
        try:
            value = self.elements[tuple]
        except KeyError:
            value = 0
        return value
