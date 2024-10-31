from common.IO import Writer, Reader

class Sphere:
    def __init__(self, coords, radius) -> None:
        if len(coords) != 3:
            raise ValueError("Coords dimension is different than 3!")
        self.coords = coords
        self.radius = radius

    def __repr__(self) -> str:
        point = ', '.join(str(coord) for coord in self.coords)
        return f"Sphere[{{{point}}},{self.radius}]"
    

class SphereContainer:
    arrSpheres = []

    def __init__(self, nMax:int = 11) -> None:
        self.max = nMax

    def __repr__(self) -> str:
        molecule = ','.join(sphere.__repr__() for sphere in self.arrSpheres).join(bracket for bracket in ['{','}'])
        return molecule 

    def getAllElements(self):
        return self.arrSpheres

    def addSphere(self, sphere: Sphere):
        self.arrSpheres.append(sphere)

    def clear(self):
        self.arrSpheres = []

    def isFull(self):
        return len(self.arrSpheres) == self.max


class SphereReader(Reader):
    def __init__(self, inputFile) -> None:
        super().__init__(inputFile)

    def readAtomElements(self):
        atom_elements = self.findLine(self.isAtomElementsHeader)
        return atom_elements[2:]

    def findLine(self, critera_func):
        line_split = self.get_line_split()
        while not critera_func(line_split):
            line_split = self.get_line_split()
        return line_split

    def isAtomElementsHeader(self, line_split):
        if len(line_split) < 5:
            return False
        return line_split[1] == 'ATOMS'


class SnapshotBuilder:
    snapshot = ''
    def addHeader(self) -> None:
        self.snapshot += "Graphics3D[{"

    def addMolecule(self, molecule: SphereContainer) -> None:
        self.snapshot += molecule.__repr__()

    def addSeparator(self, separator=',') -> None:
        self.snapshot += separator

    def addClosingBrackets(self) -> None:
        self.snapshot += "}]"

    def build(self) -> str:
        return self.snapshot


class LineParser:
    component_map = {}
    def __init__(self, atom_elements_header: list) -> None:
        for i, elem in enumerate(atom_elements_header):
            self.component_map.update({elem : i})

    def getAtomCoords(self, line):
        return [line[self.component_map[key]] for key in ('xu', 'yu', 'zu')]


class LAMMPSParser:
    def parseFile(self, fileLocation: str):
        self.setup(fileLocation)
        self.process()
        self.finalize()

    def setup(self, fileLocation: str) -> None:
        self.reader = SphereReader(fileLocation)
        self.reader.open()
        components = self.reader.readAtomElements()
        self.line = LineParser(components)
        self.container = SphereContainer()
        self.builder = SnapshotBuilder()

    def process(self) -> None:
        self.builder.addHeader()

        while not self.container.isFull():
            self.readAddSphere()

        self.builder.addMolecule(self.container)
        self.container.clear()
        i = 0
        while True:
            try:
                self.readAddSphere()
            except:
                break
            
            if not self.container.isFull():
                continue
            i += 1
            print("Molecule: ", i)
            self.builder.addSeparator()
            self.builder.addMolecule(self.container)
            self.container.clear()

        self.builder.addClosingBrackets()

    def finalize(self) -> None:
        self.reader.close()

    def readAddSphere(self) -> Sphere:
        atom = self.reader.get_line_split()
        coords = self.line.getAtomCoords(atom)
        self.container.addSphere(Sphere(coords, 0.5))

    def printSnapshot(self, targetLocation: str):
        printer = Writer(targetLocation)
        printer.open()
        snapshot = self.builder.build()
        printer.write(snapshot)
        printer.close()



if __name__ == '__main__':
    fileLocation = 'C:/Users/Szymek/Desktop/middle_snapshot_4000000.lammpstrj'
    targetLocation = 'C:/Users/Szymek/Desktop/parsed_LAMMPS.txt'

    parser = LAMMPSParser()
    parser.parseFile(fileLocation)        
    parser.printSnapshot(targetLocation)
