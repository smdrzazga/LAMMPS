from common.PARSERS import *
from common.IO import LAMMPSReader


class BananaAtom(Element):
    def __init__(self, ID, molID, type, coords) -> None:
        if len(coords) != 3:
            raise ValueError("Coords dimension is different than 3!")
        self.coords = np.array(coords)
        self.ID = int(ID)
        self.molID = int(molID)
        self.type = type

    def __repr__(self) -> str:
        return f"{self.ID} {self.molID} {self.type} {self.coords[0]} {self.coords[1]} {self.coords[2]}"
    
    

class BananaContainer(Container):
    arrElems: list[BananaAtom] = []

    def __repr__(self) -> str:
        molecule = '\n'.join(element.__repr__() for element in self.arrElems)
        return molecule 

    def changeElemTypes(self, newType):
        for i in range(len(self.arrElems)):
            self.arrElems[i].type = newType

    def updateIDs(self, offset: int):
        for i in range(len(self.arrElems)):
            self.arrElems[i].ID += offset
        
    def updateMolIDs(self, offset: int):
        for i in range(len(self.arrElems)):
            self.arrElems[i].molID += offset

    def updateTypes(self, newType: int):
        for i in range(len(self.arrElems)):
            self.arrElems[i].type = newType
        


class BananaBuilder(Builder):
    def __init__(self, sourceLocation: str):
        self.scanner = LAMMPSOutputSystemPropertiesScanner(sourceLocation)
        self.numAtoms = self.scanner.getTotalNumAtoms()
        self.box = self.scanner.getSimulationBox()

    def addHeader(self) -> None:
        header = f"""\
LAMMPS Description
 
{self.numAtoms} atoms
0 bonds
0 angles
0 dihedrals
0 impropers

1 atom types

{self.box.min[0]} {self.box.max[0]} xlo xhi
{self.box.min[1]} {self.box.max[1]} ylo yhi
{self.box.min[2]} {self.box.max[2]} zlo zhi 

Masses

1 1

Atoms 

"""
        return header
      
    def addSeparator(self) -> None:
        return "\n"

    def addClosingSymbols(self) -> None:
        return ""