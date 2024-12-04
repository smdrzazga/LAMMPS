from common.IO import Writer
from common.IO import LAMMPSReader
from common.CONTAINERS import SimulationBox
from common.WRAPPERS import *
from config import ProcessingParameters
import numpy as np
import os



class Sphere(Element):
    def __init__(self, coords, radius) -> None:
        if len(coords) != 3:
            raise ValueError("Coords dimension is different than 3!")
        self.coords = np.array(coords)
        self.radius = radius

    def __repr__(self) -> str:
        point = ', '.join(str(coord) for coord in self.coords)
        return f"Sphere[{{{point}}},{self.radius}]"


class SphereContainer(Container):
    def __repr__(self) -> str:
        molecule = ','.join(sphere.__repr__() for sphere in self.arrElems).join(bracket for bracket in ['{','}'])
        return molecule 


class SnapshotBuilder(Builder):
    def addHeader(self) -> None:
        return "Graphics3D[{"

    def addMolecule(self, molecule: Container) -> None:
        return molecule.__repr__()
        
    def addSeparator(self) -> None:
        return ','

    def addClosingSymbols(self) -> None:
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
        self.reader = ElementReader(sourceLocation)
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
            self.printer.write(self.builder.addSeparator(','))
            self.printer.write(self.builder.addMolecule(self.container))
            self.container.clear()

        self.printer.write(self.builder.addClosingSymbols())

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
        self.container.addElem(Sphere(coords, 0.5))

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
    username = os.getlogin()
    sourceLocation = 'C:/Users/' + username + '/Desktop/middle_snapshot_4000000.lammpstrj'
    targetLocation = 'C:/Users/' + username + '/Desktop/11_11p0.nb'
    boxLocation    = 'C:/Users/' + username + '/Desktop/11_11p0.ramsnap'

    parser = LAMMPSToRampackParser()
    parser.parseFile(sourceLocation, targetLocation)        
    parser.printSimulationBox(sourceLocation, boxLocation)
