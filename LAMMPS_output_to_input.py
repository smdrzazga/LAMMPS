from common.IO import Writer, Reader
from common.IO import LAMMPSReader
from common.CONTAINERS import SimulationBox
from common.WRAPPERS import *
from config import ProcessingParameters
import numpy as np
import os
from LAMMPS_to_Mathematica import *

class Atom(Element):
    def __init__(self, ID, molID, type, coords) -> None:
        if len(coords) != 3:
            raise ValueError("Coords dimension is different than 3!")
        self.coords = np.array(coords)
        self.ID = ID
        self.molID = molID
        self.type = type

    def __repr__(self) -> str:
        return f"{self.ID} {self.molID} {self.type} {self.coords[0]} {self.coords[1]} {self.coords[2]}"
    

class BananaContainer(Container):
    arrElems: list[Atom] = []

    def __repr__(self) -> str:
        molecule = '\n'.join(element.__repr__() for element in self.arrElems)
        return molecule 

    def changeElemTypes(self, newType):
        for i in range(len(self.arrElems)):
            self.arrElems[i].type = newType


class SnapshotBuilder:
    def addHeader(self) -> None:
        return "UPDATE ME PLEASE"
      
    def addSeparator(self) -> None:
        return " "

    def addClosingSymbols(self) -> None:
        return ""




if __name__ == '__main__':
    username = os.getlogin()
    sourceLocation = 'C:/Users/' + username + '/Desktop/middle_snapshot_4000000.lammpstrj'
    targetLocation = 'C:/Users/' + username + '/Desktop/11_11p0.nb'
    boxLocation    = 'C:/Users/' + username + '/Desktop/11_11p0.ramsnap'

    parser = LAMMPSToRampackParser()
    parser.parseFile(sourceLocation, targetLocation)        
    parser.printSimulationBox(sourceLocation, boxLocation)
