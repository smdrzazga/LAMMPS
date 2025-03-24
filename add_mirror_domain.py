from common.PARSERS import *
from common.BANANAS import *
from LAMMPS_output_to_input import *
import os


class MirrorBuilder(BananaBuilder):
    def __init__(self, sourceLocation):
        super().__init__(sourceLocation)
        self.numAtoms = 2*int(self.numAtoms)
        self.box.max[0] = 2.0*self.box.max[0]

    def addHeader(self) -> None:
        header = f"""\
LAMMPS Description
 
{self.numAtoms} atoms
0 bonds
0 angles
0 dihedrals
0 impropers

2 atom types

{self.box.min[0]} {self.box.max[0]} xlo xhi
{self.box.min[1]} {self.box.max[1]} ylo yhi
{self.box.min[2]} {self.box.max[2]} zlo zhi 

Masses

1 1
2 1

Atoms 
"""
        return header


class MirrorCreator(LAMMPSParser):
    mirrorParams = [120, 'x']
    moveParams = [0, 'z']

    def __init__(self, sourceLocation: str) -> None:
        self.elementClass = BananaAtom  
        self.container = BananaContainer()
        self.builder = MirrorBuilder(sourceLocation)
        self.numAtoms = int(self.builder.numAtoms)
        self.numMols  = self.numAtoms // 11


    def printSnapshot(self) -> None:
        self.printer.write(self.builder.addHeader())

        line = " "
        while True:
            try:
                line = self.readLine()
                atomParams = self.parseLine(line)
                self.addElementToContainer(*atomParams)
            except:
                break
            
            if not self.container.isFull():
                continue

            self._wrapAndPrintMolecule()
            self._wrapAndPrintMirrorMolecule()
            self.container.clear()
        self.printer.write(self.builder.addClosingSymbols())


    def parseLine(self, line) -> list[str]:
        keys = ('id', 'mol', 'element')
        atomLine = self.lineParser.getComponentByKeys(line, keys)
        atomLine.append(self.lineParser.getAtomCoords(line))
        return atomLine
        

    def _wrapAndPrintMolecule(self) -> None:
        self.container.wrapToBox(self.simBox.get_all_side_lengths())
        self.printer.write(self.builder.addSeparator())
        self.printer.write(self.builder.addMolecule(self.container))

    def _wrapAndPrintMirrorMolecule(self) -> None:
        self.container.updateIDs(self.numAtoms)
        self.container.updateMolIDs(self.numMols)
        self.container.updateTypes(2)
        self.container.mirror(*self.mirrorParams)
        self.container.move(*self.moveParams)
        self.printer.write(self.builder.addSeparator())
        self.printer.write(self.builder.addMolecule(self.container))


if __name__ == '__main__':
    username = os.getlogin()
    sourceLocation = 'C:/Users/' + username + '/Desktop/sandbox/middle_snapshot_0.lammpstrj'
    targetLocation = 'C:/Users/' + username + '/Desktop/mirror_tworzenie_atomow_z_pliku.txt'

    mirrorCreator = MirrorCreator(sourceLocation)
    mirrorCreator.processFile(sourceLocation, targetLocation)
    