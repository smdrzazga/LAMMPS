from common.IO import *
from common.PARSERS import *
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


class SphereBuilder(Builder):
    def addHeader(self) -> None:
        return "Graphics3D[{"

    def addMolecule(self, molecule: Container) -> None:
        return molecule.__repr__()
        
    def addSeparator(self) -> None:
        return ','

    def addClosingSymbols(self) -> None:
        return "}]"


class LAMMPSToMathematica(LAMMPSParser):
    def __init__(self) -> None: 
        self.elementClass = Sphere
        self.container = SphereContainer()
        self.builder = SphereBuilder()

    def parseLine(self, line) -> list[str]:
        line = [self.lineParser.getAtomCoords(line), 0.5]
        return line
    
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

    parser = LAMMPSToMathematica()
    parser.processFile(sourceLocation, targetLocation)        
    parser.printSimulationBox(sourceLocation, boxLocation)
