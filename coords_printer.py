from banana_lib import *
import numpy as np
from numpy.typing import NDArray
import os
import copy


class MoleculeParams:
    ATOMS_IN_MOL: int = 11
    CHI_DEG: float   = 110
    _ALPHA_RAD: float = (180 - CHI_DEG) * np.pi / 180
    R: float = 1 / (2 * np.sin(_ALPHA_RAD / 2 / ATOMS_IN_MOL))


def createMolecule(startAtomId: int = 1):
    molecule = Molecule(1, MoleculeParams.ATOMS_IN_MOL)

    for i in range(MoleculeParams.ATOMS_IN_MOL):
        _theta = (i - MoleculeParams.ATOMS_IN_MOL // 2) * MoleculeParams._ALPHA_RAD / (MoleculeParams.ATOMS_IN_MOL)
        coords = [MoleculeParams.R * np.cos(_theta), MoleculeParams.R * np.sin(_theta), 1.0]

        atom = Atom(startAtomId + i, coords)
        molecule.comp.append(atom)

    return molecule



class ConfigurationParameters:
    initialBoxSize: NDArray[np.float64] = np.array([60, 60, 40], dtype=np.float64)
    moleculeGrid: NDArray[np.int64]     = np.array([50, 50, 4], dtype=np.int64)
    MOL_TOTAL: int                      = np.prod(moleculeGrid)
    packingFractionEnd: float           = 0.318
    isWCA: bool                         = False
    N_WALLS: int                        = 0
    N_MOLS_PER_PERIOD: int              = 4

    offsetEmpty: NDArray[np.float64]    = np.array([8, 0, 0], dtype=np.float64)
    offsetAdd: NDArray[np.float64]      = np.array([-4, 1, 0], dtype=np.float64)
    targetFilepath: str                 = "C:/Users/" + os.getlogin() + "/Desktop/tworzenie_atomow_z_pliku.txt"



class NTBConfigurationPrinter:
    _WCAvolumeFactor = 1.015721
    _atomVolume = Atom.volume
    
    def __init__(self) -> None:
        self.volumeFactor: float = self._WCAvolumeFactor if ConfigurationParameters.isWCA else 1.0
        self.molecule: Molecule  = createMolecule()
        self.distances: NDArray  = (self.getFinalBoxSize() - ConfigurationParameters.offsetEmpty) / ConfigurationParameters.moleculeGrid


    def getFinalBoxSize(self) -> NDArray:
        _volAtoms    = ConfigurationParameters.MOL_TOTAL * self.molecule.atoms * self._atomVolume * self.volumeFactor
        _volBoxStart = np.prod(ConfigurationParameters.initialBoxSize)
        _volBoxEnd   = _volAtoms / ConfigurationParameters.packingFractionEnd

        _scale = pow(_volBoxEnd / _volBoxStart, 1 / (3 - ConfigurationParameters.N_WALLS))
        _scaleVector = [1.0, _scale, _scale] if ConfigurationParameters.N_WALLS == 1 else ([1.0, _scale, 1.0] if ConfigurationParameters.N_WALLS == 2 else [1.0, 1.0, 1.0])

        return np.array(ConfigurationParameters.initialBoxSize, dtype=np.float64) * np.array(_scaleVector, dtype=np.float64)


    def printLayer(self, layerID: int) -> None:
        _moleculeCopy = self._getRotatedMolecule(layerID)
        self._updateMoleculeProperties(_moleculeCopy)
        self._printReplicatedMolecules(_moleculeCopy, layerID)


    def _getRotatedMolecule(self, layerID) -> Molecule:
        _moleculeCopy = copy.deepcopy(self.molecule)
        position = _moleculeCopy.center_of_mass() 
        _moleculeCopy.shift(0, 0, 0)
        _moleculeCopy.rotate_x(70)
        _moleculeCopy.rotate_z(layerID * 90)
        _moleculeCopy.shift(*position)

        return _moleculeCopy


    def _updateMoleculeProperties(self, molecule: Molecule) -> None: ...


    def _printReplicatedMolecules(self, molecule: Molecule, layerID: int) -> None:
        id = ConfigurationParameters.moleculeGrid[0] * ConfigurationParameters.moleculeGrid[1] * molecule.atoms * layerID
        with open(ConfigurationParameters.targetFilepath, "a") as f:
            for i in range(ConfigurationParameters.moleculeGrid[0]):
                for j in range(ConfigurationParameters.moleculeGrid[1]):
                    for a in range(molecule.atoms):
                        id += 1
                        mol_id = (id - 1) // molecule.atoms + 1
                        itype = 1
                        x, y, z = molecule.comp[a].position + self.distances * np.array([i, j, layerID]) + ConfigurationParameters.offsetAdd
                        
                        print("%d %d %d %.5f %.2f %.5f" % (id, mol_id, itype, x, y, z), file=f)


    def printHeader(self) -> None:
        write_heading(ConfigurationParameters.targetFilepath, ConfigurationParameters.MOL_TOTAL * self.molecule.atoms, Vector(*self.getFinalBoxSize()))


    def printConfiguration(self) -> None:
        self.printHeader()
        for i in range(ConfigurationParameters.moleculeGrid[2]):
            self.printLayer(i)

    

class DipoleConfigurationPrinter(NTBConfigurationPrinter):

    def _getAtomDipoles(self, molecule: Molecule) -> NDArray[np.float64]:
        circleCenter = self._getCircleCenter(molecule)
        mainAx = self._getMainAxis(molecule)
        dipoleMoments = np.zeros((molecule.atoms, 3))

        for i, atom in enumerate(molecule.comp):
            radiusVect = (circleCenter - atom.position) / np.linalg.norm(circleCenter - atom.position)
            perpendicularAx = np.cross(mainAx, radiusVect) / np.linalg.norm(np.cross(mainAx, radiusVect))
            dipoleMoments[i] = np.cross(radiusVect, perpendicularAx)

        return dipoleMoments


    def _getCircleCenter(self, molecule: Molecule) -> NDArray[np.float64]:
        atomIndices = [1, molecule.atoms//2, -1]          # assuming N odd for simplicity
        startAtom, middleAtom, lastAtom = [molecule.comp[i].position for i in atomIndices]

        polarizationAx = (lastAtom + startAtom) / 2 - middleAtom
        polarizationAx /= np.linalg.norm(polarizationAx)

        return middleAtom + MoleculeParams.R * polarizationAx


    def _getMainAxis(self, molecule: Molecule) -> NDArray[np.float64]:
        atomIndices = [1, -1]
        startAtom, lastAtom = [molecule.comp[i].position for i in atomIndices]
        mainAx = lastAtom - startAtom
        mainAx /= np.linalg.norm(mainAx)

        return mainAx
    

    def _updateMoleculeProperties(self, molecule: Molecule) -> None:
        dipoles = self._getAtomDipoles(molecule)
        for i in range(molecule.atoms):
            molecule.comp[i].mu = dipoles[i]
            

    def _printReplicatedMolecules(self, molecule: Molecule, layerID: int) -> None:
        id = ConfigurationParameters.moleculeGrid[0] * ConfigurationParameters.moleculeGrid[1] * molecule.atoms * layerID
        with open(ConfigurationParameters.targetFilepath, "a") as f:
            for i in range(ConfigurationParameters.moleculeGrid[0]):
                for j in range(ConfigurationParameters.moleculeGrid[1]):
                    for a in range(molecule.atoms):
                        id += 1
                        mol_id = (id - 1) // molecule.atoms + 1
                        itype = 1
                        q = 0.0
                        diam = 1.0
                        rho = 1.0
                        x, y, z = molecule.comp[a].position + self.distances * np.array([i, j, layerID]) + ConfigurationParameters.offsetAdd
                        mux, muy, muz = molecule.comp[a].mu
                        #id type x y z mol q  mux muy muz  r  rho 
                        print(f"{id} {itype} {x:.5f} {y:.5f} {z:.5f} {mol_id} {q:.3f} {mux:.5f} {muy:.5f} {muz:.5f} {diam:.3f} {rho:.3f}", file=f)



if __name__ == "__main__":
    # DipoleConfigurationPrinter().printConfiguration()
    NTBConfigurationPrinter().printConfiguration()

    # molecule = createMolecule()
    # # molecule.rotate_z(180)
    # molecule.rotate_x(90)
    # molecule.rotate_y(-90)
    # molecule.rotate_x(-90)
    # molecule.shift(0, 0, 0)
    # for atom in molecule.comp:
    #     pos = atom.position
    #     pos = ", ".join([str(i) for i in pos]).join(["(", ")"])
    #     print(pos, end=",\n")

