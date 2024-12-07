from common.PARSERS import *
from common.BANANAS import *
import numpy as np
import os



class LAMMPSoutputToInput(LAMMPSParser):
    def __init__(self, sourceLocataion: str) -> None: 
        self.elementClass = BananaAtom  
        self.container = BananaContainer()
        self.builder = BananaBuilder(sourceLocataion)

    def parseLine(self, line) -> list[str]:
        keys = ('id', 'mol', 'element')
        atomLine = self.lineParser.getComponentByKeys(line, keys)
        atomLine.append(self.lineParser.getAtomCoords(line))
        return atomLine
        



if __name__ == '__main__':
    username = os.getlogin()
    sourceLocation = 'C:/Users/' + username + '/Desktop/middle_snapshot_11000000.lammpstrj'
    targetLocation = 'C:/Users/' + username + '/Desktop/test_LAMMPS_output_to_input.txt'

    parser = LAMMPSoutputToInput(sourceLocation)
    parser.processFile(sourceLocation, targetLocation)        

