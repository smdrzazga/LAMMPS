from common.IO import *
from config import ProcessingParameters
import numpy as np
from collections.abc import Iterable

class Element:
    def updatePositionAll(self, newPositions: list[float, float, float]):
        self.coords = newPositions

    def updatePosition(self, newPosition: float, axis):
        self.coords[axis] = newPosition


class Container:
    arrElems: list[Element] = []

    def __init__(self, nMax:int = 11) -> None:
        self.max = nMax

    def __repr__(self) -> str:
        raise NotImplementedError("Not implemented in the base class!")

    def getAllElements(self):
        return self.arrElems

    def addElem(self, element: Element):
        self.arrElems.append(element)

    def clear(self):
        self.arrElems = []

    def isFull(self):
        return len(self.arrElems) == self.max
    
    def mirror(self, wallPosition: float, axis: str) -> None:
        ax = {'x': 0, 'y': 1, 'z': 2}
        for i, elem in enumerate(self.arrElems):
            newPosition = 2*wallPosition - elem.coords[ax[axis]]
            self.arrElems[i].updatePosition(newPosition, ax[axis])

    def move(self, distance: float, axis: str) -> None:
        ax = {'x': 0, 'y': 1, 'z': 2}
        for i, elem in enumerate(self.arrElems):
            newPosition = elem.coords[ax[axis]] + distance
            self.arrElems[i].updatePosition(newPosition, ax[axis])

    def wrapToBox(self, boundaries: list[float, float, float]) -> None:
        atomCoords = self._getAtomCoords()
        comCoords = self._getCOMwrapped(boundaries)
        for i, atom in enumerate(self.arrElems):
            x = atomCoords[i] % boundaries
            atom.updatePositionAll(atomCoords[i] % boundaries)
        self._glueTornMolecule(comCoords, boundaries)
    
    def _getAtomCoords(self) -> list[list]:
        return np.array([element.coords for element in self.arrElems], dtype = np.float32)
    
    def _getCOMwrapped(self, boundaries: list[float, float, float]) -> list[float, float, float]:
        comCoords = np.average(self._getAtomCoords(), axis=0).flatten()
        return comCoords % boundaries
    
    def _glueTornMolecule(self, COM: list[float, float, float], boundaries: list[float, float, float]) -> None:
        _maxDist = self._maxDistToBeTorn()
        for atom in self.arrElems:
            for axis in range(3):
                dist = COM[axis] - atom.coords[axis]
                if abs(dist) < _maxDist:
                    continue    
                newPosition = atom.coords[axis] + np.sign(dist)*boundaries[axis]
                atom.updatePosition(newPosition, axis)

    def _maxDistToBeTorn(self):
        return self.max * 2


class Builder:
    def addHeader(self) -> None:
        return ''

    def addMolecule(self, molecule: Container) -> None:
        return molecule.__repr__()
        
    def addSeparator(self) -> None:
        return ''

    def addClosingSymbols(self) -> None:
        return ''
    
    
class ElementReader(Reader):
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
    

class LineParser:
    component_map = {}
    def __init__(self, atom_elements_header: list) -> None:
        for i, elem in enumerate(atom_elements_header):
            self.component_map.update({elem : i})

    def getComponentByKeys(self, line, keys: Iterable) -> list:
        return [line[self.component_map[key]] for key in keys]

    def getAtomCoords(self, line) -> list[str, str, str]:
        return self.getComponentByKeys(line, ('xu', 'yu', 'zu'))
    

class LAMMPSParser:
    def __init__(self) -> None: 
        self.elementClass = Element
        self.container = Container()
        self.builder = Builder()


    def processFile(self, sourceLocation: str, targetLocation: str) -> None:
        self.setup(sourceLocation, targetLocation)
        self.printSnapshot()
        self.finalize()


    def setup(self, sourceLocation: str, targetLocation: str) -> None:
        self.simBox = self.getSimulationBox(sourceLocation)
        self.openFiles(sourceLocation, targetLocation)
        components = self.reader.readAtomElements()
        self.lineParser = LineParser(components)


    def openFiles(self, sourceLocation: str, targetLocation: str) -> None:
        self.reader = ElementReader(sourceLocation)
        self.printer = Writer(targetLocation)
        self.reader.open()
        self.printer.open()


    def printSnapshot(self) -> None:
        self.printer.write(self.builder.addHeader())

        while not self.container.isFull():
            line = self.readLine()
            atomParams = self.parseLine(line)
            self.addElementToContainer(*atomParams)

        self.container.wrapToBox(self.simBox.get_all_side_lengths())
        self.printer.write(self.builder.addMolecule(self.container))
        self.container.clear()

        while True:
            try:
                line = self.readLine()
                atomParams = self.parseLine(line)
                self.addElementToContainer(*atomParams)
            except:
                break
            
            if not self.container.isFull():
                continue

            self.container.wrapToBox(self.simBox.get_all_side_lengths())
            self.printer.write(self.builder.addSeparator())
            self.printer.write(self.builder.addMolecule(self.container))
            self.container.clear()
        self.printer.write(self.builder.addClosingSymbols())


    def finalize(self) -> None:
        self.reader.close()
        self.printer.close()


    def getSimulationBox(self, sourceLocation) -> SimulationBox:
        reader = LAMMPSReader(ProcessingParameters(INPUT_FILE = sourceLocation, NP=1))
        reader.open(0)
        boundaries = reader.read_boundaries()
        reader.close()
        return SimulationBox(boundaries)

    def readLine(self) -> list[str]:
        return self.reader.get_line_split()

    def parseLine(self, line) -> list[str]:
        raise NotImplementedError("Function is virtual in the base class.")

    def addElementToContainer(self, *args) -> None:
        self.container.addElem(self.elementClass(*args))


class LAMMPSOutputSystemPropertiesScanner:
    def __init__(self, sourceLocation):
        self.sourceLocation = sourceLocation

    # def getTotalNumAtoms(self) -> str:
    #     with open(self.sourceLocation, 'r') as f:
    #         file = f.readlines()
    #     lastLine = file[-1]
    #     numAtoms = lastLine.strip().split()[0]
    #     return numAtoms

    def getTotalNumAtoms(self) -> str:
        with open(self.sourceLocation, 'r') as f:
            for line in f:
                if "ITEM: NUMBER " in line: 
                    numAtoms = f.readline()
                    return int(numAtoms)

        raise EOFError("The total number of atoms not found in the file!")


    def getSimulationBox(self) -> SimulationBox:
        reader = LAMMPSReader(ProcessingParameters(INPUT_FILE = self.sourceLocation, NP=1))
        reader.open(0)
        boundaries = reader.read_boundaries()
        reader.close()
        return SimulationBox(boundaries)