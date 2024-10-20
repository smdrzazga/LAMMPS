from common.IO import Writer, Reader
from config import *

class Sphere:
    def __init__(self, coords, radius) -> None:
        if len(coords) != 3:
            raise ValueError("Coords dimension is different than 3!")
        self.coords = coords
        self.radius = radius

    def __repr__(self) -> str:
        point = ', '.join(str(coord) for coord in self.coords)
        return f"Sphere[{{{point}}},{self.radius}]"

class SpherePrinter(Writer):
    def __init__(self, targetLocation) -> None:
        super().__init__(targetLocation)

    def printMolecule(self, arrSpheres: list[Sphere]):
        self.write('{', end='')
        for sphere in arrSpheres:
            self.write(sphere)
        self.write('}')

class SphereReader(Reader):
    def __init__(self, inputFile) -> None:
        super().__init__(inputFile)

    def read_atom_elements(self):
        atom_elements = self.find_line(self.is_atom_elements_header)
        return atom_elements[2:]

    def find_line(self, critera_func):
        line_split = self.get_line_split()
        while not critera_func(line_split):
            line_split = self.get_line_split()
        return self.get_line_split()

    def is_atom_elements_header(self, line_split):
        if len(line_split) < 5:
            return False
        return line_split[1] == 'ATOMS'


if __name__ == '__main__':
    fileLocation = 'C:/Users/Szymek/Desktop/middle_snapshot_4000000.lammpstrj'
    targetLocation = 'C:/Users/Szymek/Desktop/test.txt'

    reader = SphereReader(fileLocation)
    printer = SpherePrinter(targetLocation)

    reader.open()
    printer.open()

    components = reader.find_line(reader.is_atom_elements_header)
    print(components)

    printer.write(Sphere([1,2,3],1), end='\n')
    printer.printMolecule([Sphere([1,2,3],1)])

    reader.close()
    printer.close()

            


    with open(targetLocation, 'r') as t:
        outputData = t.read()
    print(outputData)
