from LAMMPS_to_Mathematica import *
import pytest
import os

class SphereTest:
    testLocation = 'C:/Users/Szymek/Desktop/test1.txt'

    def __init__(self) -> None:
        self.testSphereRepr()
        self.testPrintShpere()


    def _createTestSphere(self):
        return Sphere([1,2,3], 1)
    
    def _printToTestFile(self, text):
        printer = SpherePrinter(self.testLocation)
        printer.open()
        printer.printMolecule(text)
        printer.close()

    def _removeTestFile(self):
        if os.path.isfile(self.testLocation):
            os.remove(self.testLocation)
        else:
            raise FileNotFoundError("Error: %s file not found" % self.testLocation)
    
    def testSphereRepr(self):
        sphere = self._createTestSphere()
        assert sphere.__repr__() == "Sphere[{1, 2, 3},1]"

    def testPrintShpere(self):
        sphere = self._createTestSphere()
        self._printToTestFile([sphere])

        with open(self.testLocation, 'r') as f:
            assert f.readline() == r"{Sphere[{1, 2, 3},1]}"
        
        self._removeTestFile()