from common.SCREEN import *
from CONTAINERS import *
import numpy as np


class NTBPhase:
    def __init__(self, periods) -> None:
        self.periods = periods
        self.C_left = 0 + 0j
        self.C_right = 0 + 0j
        self.C = 0 + 0j

    def clear(self):
        self.C_left = 0 + 0j
        self.C_right = 0 + 0j
        self.C = 0 + 0j

    def update(self, molecule: Banana, box: SimulationBox) -> None:
        center = molecule.center_of_mass()
        polarization = molecule.polarization()[1]

        self.C += polarization * np.exp(2j*self.periods*np.pi * center[2] / box.get_side_length(2))
        if center[0] < box.get_side_length(0)//2:
            self.C_left += polarization * np.exp(2j*self.periods*np.pi * center[2] / box.get_side_length(2))
        else:
            self.C_right += polarization * np.exp(2j*self.periods*np.pi * center[2] / box.get_side_length(2))

    def compute_com_drift(self, box: SimulationBox):
        drift = box.get_side_length(2) / self.periods * np.angle(self.C) / (2*np.pi)
        return drift


class SmecticParameter():
    def __init__(self, periods: int) -> None:
        self.parameter = 0 + 0j
        self.count = 0
        self.periods = periods

    def __repr__(self) -> str:
       return f"{np.abs(self.parameter):.3f} | {self.count}"

    # only for vertical slices in y-z plane!
    def add_atom(self, center_coords: list, box: SimulationBox) -> None:
        self.parameter += np.exp(self.periods * 2*np.pi*1j * center_coords[-1] / box.z)
        self.count += 1

    def normalize(self) -> None:
        if self.count != 0:
            self.parameter /= self.count
    
    def read_screen(self, screen: Screen, start: int, end: int) -> None:
        HARD_CODED_LIMIT = 0

        for x in range(start, end):
            for y in range(HARD_CODED_LIMIT, screen.y):
                self.parameter += np.exp(self.periods * 2*np.pi*1j * y / screen.y) * screen.screen[y][x].colour()
                self.count += screen.screen[y][x].colour()
