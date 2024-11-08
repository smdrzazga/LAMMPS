from common.IO import Writer, Reader
from common.IO import LAMMPSReader
from common.CONTAINERS import SimulationBox
from config import ProcessingParameters
import numpy as np

class Sphere:
    def __init__(self, coords, radius) -> None:
        if len(coords) != 3:
            raise ValueError("Coords dimension is different than 3!")
        self.coords = np.array(coords)
        self.radius = radius

    def __repr__(self) -> str:
        point = ', '.join(str(coord) for coord in self.coords)
        return f"Sphere[{{{point}}},{self.radius}]"
    
    def updatePositionAll(self, newPositions: list[float, float, float]):
        self.coords = newPositions

    def updatePosition(self, newPosition: float, axis):
        self.coords[axis] = newPosition
    

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

    def wrapToBox(self, boundaries: list[float, float, float]) -> None:
        atomCoords = self._getAtomCoords()
        comCoords = self._getCOMwrapped(boundaries)
        for i, atom in enumerate(self.arrSpheres):
            x = atomCoords[i] % boundaries
            atom.updatePositionAll(atomCoords[i] % boundaries)
        self._glueTornMolecule(comCoords, boundaries)
    
    def _getAtomCoords(self) -> list[list]:
        return np.array([sphere.coords for sphere in self.arrSpheres], dtype = np.float32)
    
    def _getCOMwrapped(self, boundaries: list[float, float, float]) -> list[float, float, float]:
        comCoords = np.average(self._getAtomCoords(), axis=0).flatten()
        return comCoords % boundaries
    
    def _glueTornMolecule(self, COM: list[float, float, float], boundaries: list[float, float, float]) -> None:
        _maxDist = self.max * 2*self.arrSpheres[0].radius
        for atom in self.arrSpheres:
            for axis in range(3):
                dist = COM[axis] - atom.coords[axis]
                if abs(dist) < _maxDist:
                    continue    
                newPosition = atom.coords[axis] + np.sign(dist)*boundaries[axis]
                atom.updatePosition(newPosition, axis)


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
    def addHeader(self) -> None:
        return "Graphics3D[{"

    def addMolecule(self, molecule: SphereContainer) -> None:
        return molecule.__repr__()
        
    def addSeparator(self, separator=',') -> None:
        return separator

    def addClosingBrackets(self) -> None:
        return "}]"


class LineParser:
    component_map = {}
    def __init__(self, atom_elements_header: list) -> None:
        for i, elem in enumerate(atom_elements_header):
            self.component_map.update({elem : i})

    def getAtomCoords(self, line):
        return [line[self.component_map[key]] for key in ('xu', 'yu', 'zu')]


class LAMMPSToRampackParser:
    def parseFile(self, sourceLocation: str, targetLocation: str):
        self.setup(sourceLocation, targetLocation)
        self.printSnapshot()
        self.finalize()

    def setup(self, sourceLocation: str, targetLocation: str) -> None:
        self.simBox = self.getSimulationBox(sourceLocation)
        self.reader = SphereReader(sourceLocation)
        self.reader.open()
        self.printer = Writer(targetLocation)
        self.printer.open()
        components = self.reader.readAtomElements()
        self.line = LineParser(components)
        self.container = SphereContainer()
        self.builder = SnapshotBuilder()

    def printSnapshot(self) -> None:
        self.printer.write(self.builder.addHeader())

        while not self.container.isFull():
            self.readAddSphere()
        self.container.wrapToBox(self.simBox.get_all_side_lengths())
        self.printer.write(self.builder.addMolecule(self.container))
        self.container.clear()
        while True:
            try:
                self.readAddSphere()
            except:
                break
            
            if not self.container.isFull():
                continue
            self.container.wrapToBox(self.simBox.get_all_side_lengths())
            self.printer.write(self.builder.addSeparator())
            self.printer.write(self.builder.addMolecule(self.container))
            self.container.clear()

        self.printer.write(self.builder.addClosingBrackets())

    def finalize(self) -> None:
        self.reader.close()
        self.printer.close()

    def getSimulationBox(self, sourceLocation):
        reader = LAMMPSReader(ProcessingParameters(INPUT_FILE = sourceLocation, NP=1))
        reader.open(0)
        boundaries = reader.read_boundaries()
        reader.close()
        return SimulationBox(boundaries)

    def readAddSphere(self) -> Sphere:
        atom = self.reader.get_line_split()
        coords = self.line.getAtomCoords(atom)
        self.container.addSphere(Sphere(coords, 0.5))

    def printSimulationBox(self, sourceLocation, targetBoxLocation) -> None:
        simBox = self.getSimulationBox(sourceLocation)
        x, y, z = simBox.get_all_side_lengths()
        
        printer = Writer(targetBoxLocation)
        printer.open()
        printer.write("4", end='\n')
        printer.write("cycles LAMMPS", end='\n')
        printer.write("step.rototranslation.rotation 0.0", end='\n')
        printer.write("step.rototranslation.translation 0.0", end='\n')
        printer.write("step.scaling.scaling 0.0", end='\n')
        printer.write(str(x) + " 0.0 0.0 0.0 " + str(y) + " 0.0 0.0 0.0 " + str(z))        
        printer.close()
        



if __name__ == '__main__':
    sourceLocation = 'C:/Users/Szymek/Desktop/middle_snapshot_4000000.lammpstrj'
    targetLocation = 'C:/Users/Szymek/Desktop/11_11p0.nb'
    boxLocation    = 'C:/Users/Szymek/Desktop/11_11p0.ramsnap'

    parser = LAMMPSToRampackParser()
    parser.parseFile(sourceLocation, targetLocation)        
    parser.printSimulationBox(sourceLocation, boxLocation)
